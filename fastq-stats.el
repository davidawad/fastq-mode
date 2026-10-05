;;; fastq-stats.el --- Read and quality statistics for FASTQ files -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Package-Requires: ((emacs "29.1"))
;; Keywords: data, tools
;; URL: https://gitlab.com/davidawad/fastq-mode

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; `fastq-stats' streams records and reports read count, total bases,
;; length distribution, GC and N content, mean quality and the share of
;; bases at Q20/Q30.  By default it reads the first `fastq-stats-sample'
;; records (a few seconds); with a prefix argument it reads the whole
;; file.  When seqkit (https://bioinf.shenwei.me/seqkit/) is installed, an
;; exact whole-file summary from `seqkit stats -a' is added asynchronously.

;;; Code:

(require 'cl-lib)
(require 'fastq-core)
(require 'fastq-stream)

(defcustom fastq-stats-sample 200000
  "Records read by `fastq-stats' without a prefix argument (nil: all)."
  :type '(choice (const :tag "Whole file" nil) integer)
  :group 'fastq)

(defcustom fastq-seqkit-program "seqkit"
  "The seqkit executable, used for exact whole-file statistics if present."
  :type 'string
  :group 'fastq)

(cl-defstruct (fastq-stats (:constructor fastq-stats--make) (:copier nil))
  "Accumulated read statistics."
  (reads 0) (bases 0) (gc 0) (n 0) (qsum 0) (q20 0) (q30 0)
  (min-len most-positive-fixnum) (max-len 0)
  (lengths (make-hash-table)) offset)

(defun fastq--skip-escape (char)
  "CHAR as a literal for `skip-chars-forward'."
  (if (memq char '(?^ ?- ?\\)) (string ?\\ char) (string char)))

(defun fastq--count-run-chars (end set)
  "Count characters in SET (a `skip-chars' set) from point to END."
  (let ((n 0) (not-set (concat "^" set)))
    (while (< (point) end)
      (skip-chars-forward not-set end)
      (let ((p (point)))
        (skip-chars-forward set end)
        (cl-incf n (- (point) p))))
    n))

(defun fastq--at-least (q offset)
  "`skip-chars' set of quality characters >= Phred Q at OFFSET."
  (concat (fastq--skip-escape (+ q offset)) "-~"))

(defun fastq--stats-add (st _h s q)
  "Add the record with sequence at S and quality at Q to ST."
  (let* ((s-end (progn (goto-char s) (line-end-position)))
         (len (- s-end s))
         (q-end (progn (goto-char q) (line-end-position)))
         (off (fastq-stats-offset st))
         (h (fastq-stats-lengths st)))
    (goto-char s)
    (cl-incf (fastq-stats-gc st) (fastq--count-run-chars s-end "GCgc"))
    (goto-char s)
    (cl-incf (fastq-stats-n st) (fastq--count-run-chars s-end "Nn"))
    (cl-incf (fastq-stats-qsum st)
             (- (apply #'+ (append (buffer-substring-no-properties q q-end) nil))
                (* off (- q-end q))))
    (goto-char q)
    (cl-incf (fastq-stats-q20 st) (fastq--count-run-chars q-end (fastq--at-least 20 off)))
    (goto-char q)
    (cl-incf (fastq-stats-q30 st) (fastq--count-run-chars q-end (fastq--at-least 30 off)))
    (cl-incf (fastq-stats-reads st))
    (cl-incf (fastq-stats-bases st) len)
    (setf (fastq-stats-min-len st) (min len (fastq-stats-min-len st))
          (fastq-stats-max-len st) (max len (fastq-stats-max-len st)))
    (puthash len (1+ (gethash len h 0)) h)))

(defun fastq--file-offset (file)
  "Phred offset of FILE from its first 1000 records."
  (if (not (eq fastq-phred-offset 'auto))
      fastq-phred-offset
    (let ((text (plist-get (fastq-stream-lines file 0 4000) :text)))
      (fastq-detect-offset (mapcar #'fastq-record-quality (fastq-parse-records text))))))

(defun fastq-compute-stats (file &optional max)
  "Statistics of the first MAX records of FILE (nil: all), a `fastq-stats'."
  (let ((st (fastq-stats--make :offset (fastq--file-offset file))))
    (fastq-stream-each-record file (lambda (h s q) (fastq--stats-add st h s q)) max)
    st))

(defun fastq--pct (part whole)
  "PART as a percentage of WHOLE."
  (if (zerop whole) 0.0 (/ (* 100.0 part) whole)))

(defun fastq--length-histogram (st)
  "Lines describing the read-length distribution of ST (top 10 lengths)."
  (let* (pairs (h (fastq-stats-lengths st)) (reads (fastq-stats-reads st)))
    (maphash (lambda (k v) (push (cons k v) pairs)) h)
    (setq pairs (sort pairs (lambda (a b) (> (cdr a) (cdr b)))))
    (setq pairs (cl-subseq pairs 0 (min 10 (length pairs))))
    (mapconcat
     (lambda (p)
       (let ((pct (fastq--pct (cdr p) reads)))
         (format "  %6d bp  %10d reads  %5.1f%%  %s" (car p) (cdr p) pct
                 (make-string (round (/ pct 2.5)) ?#))))
     pairs "\n")))

(defun fastq-format-stats (st)
  "Human-readable report for ST."
  (let ((reads (fastq-stats-reads st)) (bases (fastq-stats-bases st)))
    (if (zerop reads)
        "No complete records.\n"
      (concat
       (format "Reads           %d\n" reads)
       (format "Bases           %d\n" bases)
       (format "Read length     min %d, mean %.1f, max %d\n"
               (fastq-stats-min-len st) (/ (float bases) reads) (fastq-stats-max-len st))
       (format "GC content      %.2f%%\n" (fastq--pct (fastq-stats-gc st) bases))
       (format "N (no call)     %.3f%%\n" (fastq--pct (fastq-stats-n st) bases))
       (format "Mean quality    Q%.1f (Phred+%d)\n"
               (/ (float (fastq-stats-qsum st)) (max 1 bases)) (fastq-stats-offset st))
       (format "Bases >= Q20    %.2f%%\n" (fastq--pct (fastq-stats-q20 st) bases))
       (format "Bases >= Q30    %.2f%%\n" (fastq--pct (fastq-stats-q30 st) bases))
       "\nRead lengths (most common)\n"
       (fastq--length-histogram st) "\n"))))

(defun fastq--seqkit-section (file buf)
  "Append an exact `seqkit stats -a' summary of FILE to BUF, asynchronously."
  (let ((exe (and fastq-seqkit-program (executable-find fastq-seqkit-program))))
    (when exe
      (with-current-buffer buf
        (let ((inhibit-read-only t))
          (goto-char (point-max))
          (insert "\nWhole file (seqkit stats -a): running...\n")))
      (make-process
       :name "fastq-seqkit" :buffer (generate-new-buffer " *fastq-seqkit*")
       :command (list exe "stats" "-a" "-T" (expand-file-name file))
       :connection-type 'pipe :noquery t
       :sentinel
       (lambda (proc _event)
         (unless (process-live-p proc)
           (let ((out (with-current-buffer (process-buffer proc) (buffer-string))))
             (kill-buffer (process-buffer proc))
             (when (buffer-live-p buf)
               (with-current-buffer buf
                 (let ((inhibit-read-only t))
                   (goto-char (point-max))
                   (when (search-backward "running...\n" nil t)
                     (replace-match ""))
                   (goto-char (point-max))
                   (insert (fastq--seqkit-table out))))))))))))

(defun fastq--seqkit-table (tsv)
  "Turn seqkit's two-line TSV output into aligned key/value lines."
  (let* ((lines (split-string tsv "\n" t))
         (keys (split-string (or (car lines) "") "\t"))
         (vals (split-string (or (cadr lines) "") "\t")))
    (if (< (length lines) 2)
        (concat "  seqkit failed: " (string-trim tsv) "\n")
      (mapconcat (lambda (kv) (format "  %-14s %s" (car kv) (cdr kv)))
                 (cl-mapcar #'cons keys vals) "\n"))))

;;;###autoload
(defun fastq-stats (file &optional all)
  "Show read and quality statistics for FASTQ FILE.
Reads the first `fastq-stats-sample' records, or every record with
prefix ALL; adds an exact whole-file summary when seqkit is installed."
  (interactive (list (read-file-name "FASTQ file: " nil nil t) current-prefix-arg))
  (let* ((max (if all nil fastq-stats-sample))
         (t0 (float-time))
         (st (fastq-compute-stats file max))
         (buf (get-buffer-create (format "*fastq stats %s*" (file-name-nondirectory file)))))
    (with-current-buffer buf
      (let ((inhibit-read-only t))
        (erase-buffer)
        (insert (format "%s\n%s\n\n" (abbreviate-file-name file)
                        (if (and max (= (fastq-stats-reads st) max))
                            (format "First %d records (C-u s for the whole file)" max)
                          "Whole file")))
        (insert (fastq-format-stats st))
        (insert (format "\n(%.1f s)\n" (- (float-time) t0))))
      (special-mode)
      (goto-char (point-min)))
    (fastq--seqkit-section file buf)
    (display-buffer buf)
    st))

(provide 'fastq-stats)
;;; fastq-stats.el ends here
