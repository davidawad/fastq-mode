;;; fastq-test-helpers.el --- Shared helpers for the fastq-mode tests -*- lexical-binding: t; -*-

;;; Commentary:

;; Fixture paths and generated temporary FASTQ files.  Every fixture is
;; synthetic (test/fixtures/make-fixtures.py); tests never read real data.

;;; Code:

(require 'ert)
(require 'cl-lib)

(defconst fastq-test-dir
  (file-name-directory (or load-file-name buffer-file-name))
  "The test directory.")

(defun fastq-test-fixture (name)
  "Absolute path of fixture NAME."
  (expand-file-name name (expand-file-name "fixtures" fastq-test-dir)))

(defun fastq-test-record (i &optional len)
  "Synthetic record number I (0-based) with LEN bases (default 50)."
  (let ((len (or len 50)))
    (format "@T:1:FC:1:1:%d:%d 1:N:0:AC\n%s\n+\n%s\n" i i
            (make-string len (aref "ACGT" (mod i 4)))
            (make-string len (+ 33 (mod i 41))))))

(defmacro fastq-test-with-file (spec &rest body)
  "Bind (VAR N &optional GZIP) to a temp FASTQ of N records; run BODY.
With GZIP non-nil the file is gzip-compressed (the test is skipped when
neither gzip nor zlib is available to compress it)."
  (declare (indent 1))
  (let ((var (nth 0 spec)) (n (nth 1 spec)) (gz (nth 2 spec)))
    `(let ((,var (make-temp-file "fastq-test-" nil (if ,gz ".fastq.gz" ".fastq"))))
       (unwind-protect
           (progn
             (let ((coding-system-for-write 'no-conversion))
               (with-temp-file ,var
                 (dotimes (i ,n) (insert (fastq-test-record i)))))
             (when ,gz (fastq-test--gzip-in-place ,var))
             ,@body)
         (delete-file ,var)))))

(defun fastq-test--gzip-in-place (file)
  "Compress FILE in place with gzip, or skip the current test."
  (let ((gzip (executable-find "gzip")))
    (unless gzip (ert-skip "gzip is needed to build a compressed fixture"))
    (let ((tmp (concat file ".tmp")))
      (with-temp-buffer
        (set-buffer-multibyte nil)
        (let ((coding-system-for-read 'no-conversion)
              (coding-system-for-write 'no-conversion))
          (unless (zerop (call-process gzip file t nil "-c" "-n"))
            (error "Gzip failed"))
          (write-region nil nil tmp nil 'silent)))
      (rename-file tmp file t))))

(defun fastq-test-faces-at (pos)
  "Faces at POS as a list."
  (let ((f (get-text-property pos 'face)))
    (if (listp f) f (list f))))

(provide 'fastq-test-helpers)
;;; fastq-test-helpers.el ends here
