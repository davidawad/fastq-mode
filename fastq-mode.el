;;; fastq-mode.el --- Major mode and streaming viewer for FASTQ sequencing reads -*- lexical-binding: t; -*-

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

;; FASTQ is the raw output of DNA sequencers: four lines per read (header,
;; bases, `+', per-base quality).  Files run to tens of gigabytes, usually
;; gzipped, so this package has two halves:
;;
;;   `fastq-mode'       major mode for FASTQ text you can hold in a buffer:
;;                      record-aware highlighting (bases by nucleotide,
;;                      qualities by Phred bin, low-quality bases marked),
;;                      eldoc for the base/quality at point, header decoding,
;;                      record navigation, validation.
;;   `fastq-view'       a windowed, read-only viewer for any FASTQ or
;;                      FASTQ.gz that never loads the file: it streams one
;;                      page of records at a time (see fastq-view.el).
;;
;; `fastq-open' picks the right one; `fastq-auto-view-mode' makes
;; `find-file' do the same for .fastq.gz and large .fastq files.
;; `fastq-stats' summarises reads (fastq-stats.el).  Nothing is sent
;; anywhere; external programs (gzip, seqkit) are optional.

;;; Code:

(require 'cl-lib)
(require 'fastq-core)
(require 'fastq-font)
(require 'fastq-stream)

(declare-function fastq-view "fastq-view" (file &optional record))
(declare-function fastq-view-noselect "fastq-view" (file &optional record))
(declare-function fastq-stats "fastq-stats" (file &optional all))

(defcustom fastq-view-threshold (* 20 1024 1024)
  "Plain FASTQ files larger than this many bytes open in `fastq-view'.
Compressed files always do."
  :type 'integer
  :group 'fastq)

(defconst fastq-file-name-regexp "\\.f\\(?:ast\\)?q\\(?:\\.gz\\)?\\'"
  "File names treated as FASTQ: .fastq, .fq, optionally gzipped.")

(defvar-local fastq-record-base 0
  "Number of records before the first one in this buffer (viewer pages).")

(defvar-local fastq-goto-record-function nil
  "Function taking a 1-based record number, for buffers that page records.")

(defvar-local fastq-source-file nil
  "The FASTQ file this buffer shows, for viewers whose buffer has no file.")

;;;; Records at point

(defun fastq-record-index-at-point ()
  "1-based number, within the file, of the record at point."
  (+ fastq-record-base 1 (/ (1- (line-number-at-pos)) fastq-lines-per-record)))

(defun fastq--record-start ()
  "Position of the header line of the record at point."
  (save-excursion
    (forward-line (- (mod (1- (line-number-at-pos)) fastq-lines-per-record)))
    (point)))

(defun fastq-record-at-point ()
  "The record at point as a `fastq-record', or nil if it is incomplete."
  (save-excursion
    (goto-char (fastq--record-start))
    (let (lines)
      (dotimes (_ fastq-lines-per-record)
        (push (buffer-substring-no-properties (point) (line-end-position)) lines)
        (forward-line 1))
      (pcase-let ((`(,q ,_ ,s ,h) lines))
        (and q (string-prefix-p "@" h)
             (fastq-record-create :header h :sequence s :quality q))))))

(defun fastq--detect-buffer-offset ()
  "Phred offset of the qualities in the first records of this buffer."
  (save-excursion
    (goto-char (point-min))
    (let (quals)
      (while (and (< (length quals) 1000) (zerop (forward-line 3)) (not (eobp)))
        (push (buffer-substring-no-properties (point) (line-end-position)) quals)
        (forward-line 1))
      (fastq-offset quals))))

;;;; Eldoc

(defun fastq--base-doc (rec i)
  "Eldoc text for base I (0-based) of record REC."
  (let* ((seq (fastq-record-sequence rec)) (qual (fastq-record-quality rec))
         (off fastq-buffer-offset))
    (if (>= i (min (length seq) (length qual)))
        (format "Read %d: %d bases, mean Q%.1f"
                (fastq-record-index-at-point) (length seq)
                (fastq-mean-quality qual off))
      (let ((q (fastq-phred (aref qual i) off)))
        (format "Read %d, base %d/%d: %c  Q%d (error %.3g%%, %s)"
                (fastq-record-index-at-point) (1+ i) (length seq) (aref seq i) q
                (* 100 (fastq-error-probability q)) (fastq-quality-bin q))))))

(defun fastq-eldoc-function (&rest _)
  "Describe the header, base or quality at point."
  (let ((rec (fastq-record-at-point))
        (col (- (point) (line-beginning-position))))
    (when rec
      (pcase (fastq-line-role (line-number-at-pos))
        ('header (format "Read %d: %s" (fastq-record-index-at-point)
                         (fastq-header-summary (fastq-record-header rec))))
        ((or 'sequence 'quality) (fastq--base-doc rec col))
        ('separator (format "Read %d: %d bases" (fastq-record-index-at-point)
                            (length (fastq-record-sequence rec))))))))

;;;; Commands

(defun fastq-next-record (&optional n)
  "Move to the header of the Nth next record."
  (interactive "p")
  (goto-char (fastq--record-start))
  (forward-line (* fastq-lines-per-record (or n 1))))

(defun fastq-previous-record (&optional n)
  "Move to the header of the Nth previous record."
  (interactive "p")
  (fastq-next-record (- (or n 1))))

(defun fastq-goto-record (n)
  "Move to record N (1-based) of the file."
  (interactive "nGo to record: ")
  (if fastq-goto-record-function
      (funcall fastq-goto-record-function n)
    (goto-char (point-min))
    (forward-line (* fastq-lines-per-record (1- (max 1 n))))))

(defun fastq-current-file ()
  "The FASTQ file shown in this buffer."
  (or fastq-source-file buffer-file-name
      (user-error "This buffer is not visiting a FASTQ file")))

(defun fastq-describe-header ()
  "Show every decoded field of the read header at point."
  (interactive)
  (let ((rec (or (fastq-record-at-point) (user-error "No complete record at point"))))
    (with-help-window "*fastq header*"
      (princ (fastq-record-header rec))
      (princ "\n\n")
      (pcase-dolist (`(,k . ,v) (fastq-decode-header (fastq-record-header rec)))
        (princ (format "%-11s %s\n" k v))))))

(defun fastq-validate ()
  "Check every record in the buffer; stop at the first problem."
  (interactive)
  (goto-char (point-min))
  (let ((n 0) problem)
    (while (and (not problem) (not (eobp)))
      (let ((lines (cl-loop repeat 4
                            collect (prog1 (buffer-substring-no-properties
                                            (point) (line-end-position))
                                      (forward-line 1)))))
        (cl-incf n)
        (setq problem (apply #'fastq-validate-record lines))))
    (if problem
        (progn (fastq-previous-record)
               (message "Record %d: %s" (+ fastq-record-base n) problem))
      (message "%d records, all valid" n))))

(defun fastq-mate (&optional record)
  "Show the paired-end mate file (R1 <-> R2) at RECORD (default: this one)."
  (interactive)
  (let* ((file (fastq-current-file))
         (mate (or (fastq-mate-file-name file)
                   (user-error "No R1/R2 pattern in %s" (file-name-nondirectory file)))))
    (unless (file-exists-p mate) (user-error "Mate file %s not found" mate))
    (let ((buf (fastq-open-noselect mate (or record (fastq-record-index-at-point))))
          (name (and (fastq-record-at-point)
                     (fastq-read-name (fastq-record-header (fastq-record-at-point))))))
      (pop-to-buffer buf)
      (let ((other (fastq-record-at-point)))
        (when (and name other
                   (not (equal name (fastq-read-name (fastq-record-header other)))))
          (message "Warning: mate read name %s differs from %s"
                   (fastq-read-name (fastq-record-header other)) name))))))

(defun fastq-stats-this-file (&optional all)
  "Summarise this buffer's FASTQ file; with prefix ALL, read every record."
  (interactive "P")
  (require 'fastq-stats)
  (fastq-stats (fastq-current-file) all))

;;;; Opening files

(defun fastq-view-file-p (file)
  "Non-nil if FILE should open in the streaming viewer."
  (or (fastq-gzip-file-p file)
      (> (or (file-attribute-size (file-attributes file)) 0) fastq-view-threshold)))

(defun fastq-open-noselect (file &optional record)
  "Return a buffer showing FASTQ FILE at RECORD, viewer or `fastq-mode'."
  (require 'fastq-view)
  (if (fastq-view-file-p file)
      (fastq-view-noselect file record)
    (let ((buf (find-file-noselect file)))
      (with-current-buffer buf
        (unless (derived-mode-p 'fastq-mode) (fastq-mode))
        (when record (fastq-goto-record record)))
      buf)))

;;;###autoload
(defun fastq-open (file &optional record)
  "Open FASTQ FILE, streaming it when it is compressed or large.
With a prefix argument, ask for the RECORD number to start at."
  (interactive (list (read-file-name "FASTQ file: " nil nil t)
                     (and current-prefix-arg (read-number "Record: " 1))))
  (switch-to-buffer (fastq-open-noselect file record)))

(defun fastq--find-file-advice (orig filename &rest args)
  "Around advice for `find-file-noselect': stream FASTQ files that need it.
ORIG is the advised function, FILENAME and ARGS its arguments."
  (let ((f (expand-file-name filename)))
    (if (and (string-match-p fastq-file-name-regexp f)
             (file-regular-p f) (file-readable-p f)
             (fastq-view-file-p f))
        (progn (require 'fastq-view) (fastq-view-noselect f))
      (apply orig filename args))))

;;;###autoload
(define-minor-mode fastq-auto-view-mode
  "Open .fastq.gz and large .fastq files in the streaming `fastq-view'.
Without it, `find-file' would decompress or read the whole file."
  :global t :group 'fastq
  (if fastq-auto-view-mode
      (advice-add 'find-file-noselect :around #'fastq--find-file-advice)
    (advice-remove 'find-file-noselect #'fastq--find-file-advice)))

;;;; Mode

(defvar-keymap fastq-mode-map
  :doc "Keymap for `fastq-mode'."
  "C-c C-n" #'fastq-next-record
  "C-c C-p" #'fastq-previous-record
  "C-c C-g" #'fastq-goto-record
  "C-c C-h" #'fastq-describe-header
  "C-c C-m" #'fastq-mate
  "C-c C-s" #'fastq-stats-this-file
  "C-c C-v" #'fastq-validate)

;;;###autoload
(define-derived-mode fastq-mode text-mode "FASTQ"
  "Major mode for FASTQ sequencing reads.
Bases are coloured by nucleotide and qualities by Phred bin; eldoc
describes the base or header at point.  For compressed or very large
files use `fastq-open', which streams them instead.

\\{fastq-mode-map}"
  (setq-local truncate-lines t)
  (setq-local fastq-buffer-offset (fastq--detect-buffer-offset))
  (fastq-font-setup)
  (add-hook 'eldoc-documentation-functions #'fastq-eldoc-function nil t)
  (eldoc-mode 1))

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.f\\(?:ast\\)?q\\'" . fastq-mode))

(provide 'fastq-mode)
;;; fastq-mode.el ends here
