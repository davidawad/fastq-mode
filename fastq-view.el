;;; fastq-view.el --- Windowed streaming viewer for huge FASTQ files -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Keywords: data, tools
;; URL: https://gitlab.com/davidawad/fastq-mode

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; `fastq-view' shows one page of `fastq-view-page-records' records from a
;; FASTQ or FASTQ.gz file of any size; the rest of the file is never read
;; into Emacs.  Paging forward streams on from the file; jumping to record
;; N in a plain file seeks from the nearest remembered byte offset (one is
;; kept every `fastq-checkpoint-interval' records, per file, for the
;; session), while a gzip file has to be decompressed up to N.
;;
;; Keys: n/p record, SPC/DEL page, g goto record, < > first/last page,
;; # count records, m mate (R1<->R2), s stats, h header, w copy record,
;; v validate page, q quit.

;;; Code:

(require 'cl-lib)
(require 'fastq-mode)

(defcustom fastq-view-page-records 200
  "Records shown per page in `fastq-view'."
  :type 'integer
  :group 'fastq)

(defvar fastq-view--checkpoints (make-hash-table :test #'equal)
  "File identity -> alist of (RECORD . BYTE) for plain files.")

(defvar-local fastq-view-total nil
  "Number of records in the file, once known.")

(defun fastq-view--file-key (file)
  "Identity of FILE for the checkpoint cache: name, size and mtime."
  (let ((a (file-attributes file)))
    (list (file-truename file) (file-attribute-size a)
          (file-attribute-modification-time a))))

(defun fastq-view--nearest-checkpoint (file record)
  "Best known (RECORD0 . BYTE) at or before 0-based RECORD of FILE."
  (let ((best '(0 . 0)))
    (dolist (cp (gethash (fastq-view--file-key file) fastq-view--checkpoints))
      (when (and (<= (car cp) record) (> (car cp) (car best)))
        (setq best cp)))
    best))

(defun fastq-view--remember (file cp)
  "Remember checkpoint CP, (RECORD . BYTE), for FILE."
  (let ((key (fastq-view--file-key file)))
    (unless (assq (car cp) (gethash key fastq-view--checkpoints))
      (push cp (gethash key fastq-view--checkpoints)))))

(defun fastq-view-fetch (file start count)
  "Fetch COUNT records of FILE from 0-based record START.
Return the `fastq-stream-lines' plist, plus :start, the record actually
shown first (moved back when START is past the end of the file)."
  (let* ((plain (eq (fastq-file-kind file) 'plain))
         (cp (if plain (fastq-view--nearest-checkpoint file start) '(0 . 0)))
         (res (fastq-stream-lines
               file (* 4 (- start (car cp))) (* 4 count)
               :start-byte (cdr cp) :first-record (car cp)
               :on-checkpoint (and plain (lambda (c) (fastq-view--remember file c))))))
    (when (plist-get res :byte)
      (fastq-view--remember file (cons start (plist-get res :byte))))
    (if (and (zerop (plist-get res :lines)) (> start 0))
        ;; past the end: show the last page instead
        (let ((total (+ (car cp) (/ (plist-get res :skipped) 4))))
          (append (fastq-view-fetch file (max 0 (- total count)) count)
                  (list :total total)))
      (append res (list :start start)))))

(defun fastq-view--header-line ()
  "Header line for the viewer buffer."
  (let* ((n (/ (count-lines (point-min) (point-max)) 4))
         (first (1+ fastq-record-base)) (last (+ fastq-record-base n)))
    (format " %s   records %s-%s of %s   Phred+%d   n/p record  SPC/DEL page  g goto  m mate  s stats  h header  q quit"
            (file-name-nondirectory fastq-source-file)
            first last
            (if fastq-view-total (number-to-string fastq-view-total) "?")
            fastq-buffer-offset)))

(defun fastq-view--load (start &optional at-end)
  "Show the page starting at 0-based record START.
Point goes to the first record, or the last one when AT-END."
  (let* ((res (fastq-view-fetch fastq-source-file (max 0 start) fastq-view-page-records))
         (inhibit-read-only t))
    (when (plist-get res :total) (setq fastq-view-total (plist-get res :total)))
    (when (and (plist-get res :eof) (not fastq-view-total))
      (setq fastq-view-total (+ (plist-get res :start) (/ (plist-get res :lines) 4))))
    (erase-buffer)
    (insert (plist-get res :text))
    (setq fastq-record-base (plist-get res :start))
    (goto-char (point-min))
    (when at-end
      (goto-char (point-max))
      (forward-line -4))
    (set-buffer-modified-p nil)
    (setq header-line-format '(:eval (fastq-view--header-line)))))

(defun fastq-view--records-on-page ()
  "Number of whole records in the buffer."
  (/ (count-lines (point-min) (point-max)) 4))

(defun fastq-view-goto-record (n)
  "Show record N (1-based), loading its page if needed."
  (interactive "nGo to record: ")
  (let ((i (1- (max 1 n))))
    (cond
     ((< i fastq-record-base)
      ;; going back: the page ends at the record
      (fastq-view--load (max 0 (- i (1- fastq-view-page-records)))))
     ((>= i (+ fastq-record-base (fastq-view--records-on-page)))
      (fastq-view--load i)))
    (goto-char (point-min))
    (forward-line (* 4 (max 0 (min (- i fastq-record-base)
                                   (1- (fastq-view--records-on-page))))))))

(defun fastq-view-next-record (&optional n)
  "Move N records forward, crossing pages."
  (interactive "p")
  (fastq-view-goto-record (+ (fastq-record-index-at-point) (or n 1))))

(defun fastq-view-previous-record (&optional n)
  "Move N records back, crossing pages."
  (interactive "p")
  (fastq-view-goto-record (max 1 (- (fastq-record-index-at-point) (or n 1)))))

(defun fastq-view-next-page ()
  "Show the next page of records."
  (interactive)
  (let ((next (+ fastq-record-base (fastq-view--records-on-page))))
    (if (and fastq-view-total (>= next fastq-view-total))
        (message "End of %s" (file-name-nondirectory fastq-source-file))
      (fastq-view--load next))))

(defun fastq-view-previous-page ()
  "Show the previous page of records."
  (interactive)
  (if (zerop fastq-record-base)
      (message "Beginning of %s" (file-name-nondirectory fastq-source-file))
    (fastq-view--load (max 0 (- fastq-record-base fastq-view-page-records)) t)))

(defun fastq-view-first-page ()
  "Show the first page."
  (interactive)
  (fastq-view--load 0))

(defun fastq-view-count ()
  "Count the records in the file (streams all of it once)."
  (interactive)
  (setq fastq-view-total (fastq-stream-count-records fastq-source-file))
  (force-mode-line-update)
  (message "%s: %d records" (file-name-nondirectory fastq-source-file) fastq-view-total))

(defun fastq-view-last-page ()
  "Show the last page (counts the records first if needed)."
  (interactive)
  (unless fastq-view-total (fastq-view-count))
  (fastq-view--load (max 0 (- fastq-view-total fastq-view-page-records)) t))

(defun fastq-view-copy-record ()
  "Copy the record at point to the kill ring."
  (interactive)
  (let ((rec (or (fastq-record-at-point) (user-error "No record at point"))))
    (kill-new (format "%s\n%s\n+\n%s\n" (fastq-record-header rec)
                      (fastq-record-sequence rec) (fastq-record-quality rec)))
    (message "Copied read %d" (fastq-record-index-at-point))))

(defun fastq-view-revert (&rest _)
  "Reload the current page from the file."
  (fastq-view--load fastq-record-base))

(defvar-keymap fastq-view-mode-map
  :doc "Keymap for `fastq-view-mode'."
  :parent special-mode-map
  "n" #'fastq-view-next-record
  "p" #'fastq-view-previous-record
  "SPC" #'fastq-view-next-page
  "DEL" #'fastq-view-previous-page
  "S-SPC" #'fastq-view-previous-page
  "g" #'fastq-view-goto-record
  "<" #'fastq-view-first-page
  ">" #'fastq-view-last-page
  "#" #'fastq-view-count
  "m" #'fastq-mate
  "s" #'fastq-stats-this-file
  "h" #'fastq-describe-header
  "v" #'fastq-validate
  "w" #'fastq-view-copy-record)

(define-derived-mode fastq-view-mode fastq-mode "FASTQ-view"
  "Read-only, paged view of a FASTQ file that is streamed, never loaded.

\\{fastq-view-mode-map}"
  (setq buffer-read-only t)
  (setq-local revert-buffer-function #'fastq-view-revert)
  (setq-local fastq-goto-record-function #'fastq-view-goto-record))

;;;###autoload
(defun fastq-view-noselect (file &optional record)
  "Return a `fastq-view' buffer for FILE showing RECORD (1-based)."
  (let* ((file (expand-file-name file))
         (name (format "*fastq %s*" (file-name-nondirectory file)))
         (buf (or (cl-find-if (lambda (b) (equal (buffer-local-value 'fastq-source-file b) file))
                              (buffer-list))
                  (generate-new-buffer name))))
    (with-current-buffer buf
      (unless (derived-mode-p 'fastq-view-mode)
        (fastq-view-mode)
        (setq fastq-source-file file)
        (fastq-view--load 0)
        (setq fastq-buffer-offset (fastq--detect-buffer-offset))
        (font-lock-flush))
      (setq default-directory (file-name-directory file))
      (when record (fastq-view-goto-record record)))
    buf))

;;;###autoload
(defun fastq-view (file &optional record)
  "Browse FASTQ FILE (plain or gzipped, any size) a page at a time.
With a prefix argument, ask for the RECORD number to start at."
  (interactive (list (read-file-name "FASTQ file: " nil nil t)
                     (and current-prefix-arg (read-number "Record: " 1))))
  (switch-to-buffer (fastq-view-noselect file record)))

(provide 'fastq-view)
;;; fastq-view.el ends here
