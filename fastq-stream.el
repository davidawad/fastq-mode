;;; fastq-stream.el --- Read FASTQ lines without loading the file -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Keywords: data, tools
;; URL: https://github.com/davidawad/fastq-mode

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; A "source" feeds bytes of a FASTQ file into a unibyte scratch buffer a
;; chunk at a time, so memory stays bounded however large the file is:
;;
;;   plain    `insert-file-contents-literally' with byte ranges (seekable)
;;   gzip     `gzip -dc' through a pipe; the process is killed as soon as
;;            enough lines have arrived (gzip cannot seek, so reaching
;;            record N costs a decompression scan of the records before it)
;;   zlib     Emacs' built-in zlib, for small .gz files on machines with
;;            no gzip executable (it must decompress the whole file)
;;
;; `fastq-stream-lines' skips and takes lines; `fastq-stream-each-record'
;; walks records for statistics.  For plain files a caller-supplied
;; checkpoint function receives (RECORD . BYTE) pairs so later jumps can
;; seek instead of scanning.

;;; Code:

(require 'cl-lib)
(require 'fastq-core)

(defcustom fastq-gzip-program "gzip"
  "Program used to decompress .gz files (run as PROGRAM -dc -- FILE).
pigz or a gzip from Git for Windows work too.  When it is not found,
small files fall back to Emacs' own zlib (see `fastq-zlib-max-size')."
  :type 'string
  :group 'fastq)

(defcustom fastq-zlib-max-size (* 64 1024 1024)
  "Largest compressed file, in bytes, that may be decompressed with zlib.
zlib inflates the whole file in memory, so it is only a fallback for
machines without `fastq-gzip-program'."
  :type 'integer
  :group 'fastq)

(defcustom fastq-chunk-size (* 4 1024 1024)
  "Bytes read from a plain file per chunk."
  :type 'integer
  :group 'fastq)

(defcustom fastq-checkpoint-interval 10000
  "Records between remembered byte offsets in plain (uncompressed) files."
  :type 'integer
  :group 'fastq)

;; Defined only in builds with zlib.
(declare-function zlib-available-p "decompress.c" ())
(declare-function zlib-decompress-region "decompress.c"
                  (start end &optional allow-partial))

(defun fastq-gzip-file-p (file)
  "Return non-nil if FILE begins with the gzip magic bytes."
  (with-temp-buffer
    (set-buffer-multibyte nil)
    (insert-file-contents-literally file nil 0 2)
    (equal (buffer-string) "\x1f\x8b")))

(defun fastq--gzip-executable ()
  "The gzip executable, or nil."
  (and fastq-gzip-program (executable-find fastq-gzip-program)))

(defun fastq--zlib-p ()
  "Non-nil if this Emacs can inflate gzip data itself."
  (and (fboundp 'zlib-available-p) (zlib-available-p)))

(defun fastq-file-kind (file)
  "How FILE is read: `plain', `gzip' or `zlib'.
Signals `fastq-error' if a .gz file cannot be decompressed here."
  (cond
   ((not (fastq-gzip-file-p file)) 'plain)
   ((fastq--gzip-executable) 'gzip)
   ((and (fastq--zlib-p)
         (<= (file-attribute-size (file-attributes file)) fastq-zlib-max-size))
    'zlib)
   (t (fastq--error
       "Cannot stream %s: install gzip (or set `fastq-gzip-program'); zlib is limited to files under %s"
       (file-name-nondirectory file) (file-size-human-readable fastq-zlib-max-size)))))

;;;; Sources

(cl-defstruct (fastq-source (:constructor fastq-source--make) (:copier nil))
  "A chunked reader filling the current buffer."
  file kind proc errbuf (pos 0) (extra 0) eof)

(defun fastq-source-open (file &optional start-byte)
  "Open FILE for reading into the current (unibyte) buffer.
START-BYTE seeks plain files; it must be 0 or nil for compressed ones."
  (let* ((kind (fastq-file-kind file))
         (src (fastq-source--make :file file :kind kind :pos (or start-byte 0))))
    (when (and (not (eq kind 'plain)) (> (fastq-source-pos src) 0))
      (fastq--error "Cannot seek in compressed file %s" file))
    (when (eq kind 'gzip)
      (let ((err (generate-new-buffer " *fastq-gunzip-stderr*")))
        (setf (fastq-source-errbuf src) err
              (fastq-source-proc src)
              (make-process
               :name "fastq-gunzip" :buffer (current-buffer)
               :command (list (fastq--gzip-executable) "-dc" "--"
                              (expand-file-name file))
               :connection-type 'pipe :coding 'binary :noquery t
               :stderr err :sentinel #'ignore))))
    src))

(defun fastq--fill-plain (src)
  "Append the next chunk of plain SRC; nil at end of file."
  (goto-char (point-max))
  (let* ((pos (fastq-source-pos src))
         (n (cadr (insert-file-contents-literally
                   (fastq-source-file src) nil pos (+ pos fastq-chunk-size)))))
    (cl-incf (fastq-source-pos src) n)
    (> n 0)))

(defun fastq--fill-zlib (src)
  "Inflate all of SRC (zlib fallback) into the buffer, once."
  (unless (> (fastq-source-pos src) 0)
    (goto-char (point-max))
    (let ((beg (point)))
      (insert-file-contents-literally (fastq-source-file src))
      (goto-char (point-max))
      (unless (zlib-decompress-region beg (point-max))
        (fastq--error "Corrupt gzip data in %s" (fastq-source-file src))))
    (setf (fastq-source-pos src) 1)
    t))

(defun fastq--fill-gzip (src)
  "Wait for more output from the gzip process of SRC; nil once it is done."
  (let ((proc (fastq-source-proc src))
        (size (buffer-size))
        (read-process-output-max (* 1024 1024)))
    (while (and (= size (buffer-size))
                (or (accept-process-output proc 0.2)
                    (process-live-p proc))))
    (when (and (= size (buffer-size)) (not (process-live-p proc)))
      (while (accept-process-output proc 0))
      (unless (or (> (buffer-size) size) (zerop (process-exit-status proc)))
        (fastq--error "Gzip failed on %s (exit %s): %s" (fastq-source-file src)
                      (process-exit-status proc)
                      (with-current-buffer (fastq-source-errbuf src)
                        (string-trim (buffer-string))))))
    (> (buffer-size) size)))

(defun fastq-source-fill (src)
  "Append more data from SRC to the current buffer; nil at end of input.
At the end, a missing final newline is supplied (once) so the last line
counts as complete."
  (unless (fastq-source-eof src)
    (or (pcase (fastq-source-kind src)
          ('plain (fastq--fill-plain src))
          ('zlib (fastq--fill-zlib src))
          ('gzip (fastq--fill-gzip src)))
        (progn
          (setf (fastq-source-eof src) t)
          (when (and (> (buffer-size) 0) (/= (char-before (point-max)) ?\n))
            (goto-char (point-max))
            (insert "\n")
            (cl-incf (fastq-source-extra src))
            t)))))

(defun fastq-source-byte (src)
  "Byte offset in the file of the current buffer's first character (plain SRC)."
  (- (+ (fastq-source-pos src) (fastq-source-extra src)) (buffer-size)))

(defun fastq-source-close (src)
  "Stop SRC's processes and drop its stderr buffer."
  (let ((p (fastq-source-proc src)) (err (fastq-source-errbuf src)))
    (when (and p (process-live-p p)) (delete-process p))
    (when (buffer-live-p err)
      (let ((ep (get-buffer-process err)))
        (when ep (delete-process ep)))
      (kill-buffer err))))

;;;; Lines

(defun fastq--forward-complete-lines (n)
  "Move over up to N complete (newline-terminated) lines; return the count."
  (let* ((short (forward-line n)) (moved (- n short)))
    (if (and (> moved 0) (not (bolp)))
        (progn (forward-line 0) (1- moved))
      moved)))

(defun fastq--skip-lines (src skip on-checkpoint first-record)
  "Discard SKIP lines from the front of the buffer, filling from SRC.
Return the number of lines actually skipped (less at end of input).
ON-CHECKPOINT, if non-nil, is called with (RECORD . BYTE) at every
`fastq-checkpoint-interval' records of a plain SRC; FIRST-RECORD is the
0-based record number of the buffer's first line."
  (let ((left skip) (base (fastq-source-byte src))
        (lines-done 0) (last-msg (float-time)))
    (while (> left 0)
      (goto-char (point-min))
      (let* ((to-cp (and on-checkpoint (eq (fastq-source-kind src) 'plain)
                         (let* ((rec (+ first-record (/ lines-done 4)))
                                (next (* fastq-checkpoint-interval
                                         (1+ (/ rec fastq-checkpoint-interval)))))
                           (- (* 4 (- next first-record)) lines-done))))
             (step (if to-cp (min left to-cp) left))
             (moved (fastq--forward-complete-lines step)))
        (cl-incf base (1- (point)))
        (delete-region (point-min) (point))
        (cl-decf left moved)
        (cl-incf lines-done moved)
        (when (and to-cp (= moved to-cp))
          (funcall on-checkpoint (cons (+ first-record (/ lines-done 4)) base)))
        (when (and (> left 0) (< moved step) (not (fastq-source-fill src)))
          (setq left 0 skip lines-done))
        (when (> (- (float-time) last-msg) 1.0)
          (setq last-msg (float-time))
          (message "fastq: scanning %s... %d records"
                   (file-name-nondirectory (fastq-source-file src))
                   (+ first-record (/ lines-done 4))))))
    (min skip lines-done)))

(cl-defun fastq-stream-lines (file skip take &key (start-byte 0) (first-record 0)
                                   on-checkpoint)
  "Return lines SKIP .. SKIP+TAKE of FILE without loading the rest of it.
START-BYTE (plain files only) is where reading begins and FIRST-RECORD
the record number found there.  The result is a plist:
  :text     the lines, newline-terminated
  :lines    how many lines :text holds (fewer than TAKE at end of file)
  :skipped  lines actually skipped (fewer than SKIP if the file is shorter)
  :byte     byte offset of :text in a plain file, else nil
  :eof      non-nil if the file ended within or before :text
ON-CHECKPOINT is passed (RECORD . BYTE) pairs; see `fastq--skip-lines'."
  (with-temp-buffer
    (set-buffer-multibyte nil)
    (let ((src (fastq-source-open file start-byte)))
      (unwind-protect
          (let* ((skipped (fastq--skip-lines src skip on-checkpoint first-record))
                 (byte (and (eq (fastq-source-kind src) 'plain)
                            (fastq-source-byte src)))
                 (got 0))
            (while (progn (goto-char (point-min))
                          (setq got (fastq--forward-complete-lines take))
                          (and (< got take) (fastq-source-fill src))))
            (list :text (buffer-substring-no-properties (point-min) (point))
                  :lines got :skipped skipped :byte byte
                  :eof (fastq-source-eof src)))
        (fastq-source-close src)))))

(defun fastq-stream-count-records (file)
  "Count the records in FILE by streaming all of it."
  (/ (plist-get (fastq-stream-lines file most-positive-fixnum 0) :skipped) 4))

;;;; Records

(defun fastq-stream-each-record (file fn &optional max)
  "Call FN on each record of FILE, at most MAX records (nil: all).
FN is called in a scratch buffer with point at the record's header and
the buffer positions (HEADER SEQ QUAL) of the starts of its header,
sequence and quality lines; each line ends at `line-end-position' from
there.  Return the number of records visited."
  (with-temp-buffer
    (set-buffer-multibyte nil)
    (let ((src (fastq-source-open file)) (n 0) (more t))
      (unwind-protect
          (while (and more (or (null max) (< n max)))
            (goto-char (point-min))
            (let ((done nil))
              (while (and (not done) (or (null max) (< n max)))
                (let ((start (point)))
                  (if (/= (fastq--forward-complete-lines 4) 4)
                      (progn (goto-char start) (setq done t))
                    (let ((end (point)))
                      (goto-char start)
                      (let* ((h (point))
                             (s (progn (forward-line 1) (point)))
                             (q (progn (forward-line 2) (point))))
                        (save-excursion (goto-char h) (funcall fn h s q)))
                      (goto-char end)
                      (cl-incf n)))))
              (delete-region (point-min) (point))
              (setq more (fastq-source-fill src))))
        (fastq-source-close src))
      n)))

(provide 'fastq-stream)
;;; fastq-stream.el ends here
