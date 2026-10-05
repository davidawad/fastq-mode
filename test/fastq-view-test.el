;;; fastq-view-test.el --- Tests for the streaming viewer -*- lexical-binding: t; -*-

;;; Commentary:

;; Paging, jumping across pages, end-of-file handling, counting, copying,
;; and the plain-file checkpoint cache.

;;; Code:

(require 'fastq-test-helpers)
(require 'fastq-view)

(defmacro fastq-view-test-with (spec &rest body)
  "Open (FILE &optional PAGE-RECORDS) in a viewer buffer and run BODY there."
  (declare (indent 1))
  `(let* ((fastq-view-page-records (or ,(nth 1 spec) 10))
          (buf (fastq-view-noselect ,(nth 0 spec))))
     (unwind-protect (with-current-buffer buf ,@body)
       (kill-buffer buf))))

(defun fastq-view-test--header-at-point ()
  "Header line of the record at point."
  (fastq-record-header (fastq-record-at-point)))

(ert-deftest fastq-view-test-first-page ()
  (fastq-test-with-file (f 95)
    (fastq-view-test-with (f 10)
      (should (eq major-mode 'fastq-view-mode))
      (should buffer-read-only)
      (should (= (count-lines (point-min) (point-max)) 40))
      (should (= (fastq-record-index-at-point) 1))
      (should (string-match-p "records 1-10 of \\?" (fastq-view--header-line))))))

(ert-deftest fastq-view-test-paging-and-records ()
  (fastq-test-with-file (f 95)
    (fastq-view-test-with (f 10)
      (fastq-view-next-page)
      (should (= fastq-record-base 10))
      (should (string-match-p ":1:1:10:" (fastq-view-test--header-at-point)))
      ;; n across the page boundary loads the next page
      (fastq-view-goto-record 20)
      (fastq-view-next-record)
      (should (= fastq-record-base 20))
      (should (= (fastq-record-index-at-point) 21))
      ;; p across the boundary goes back
      (fastq-view-previous-record)
      (should (= (fastq-record-index-at-point) 20))
      (should (string-match-p ":1:1:19:" (fastq-view-test--header-at-point)))
      (fastq-view-previous-page)
      (should (= fastq-record-base 0)))))

(ert-deftest fastq-view-test-end-of-file ()
  (fastq-test-with-file (f 95)
    (fastq-view-test-with (f 10)
      (fastq-view-goto-record 500)
      ;; past the end shows the last page and learns the total
      (should (= fastq-view-total 95))
      (should (= fastq-record-base 85))
      (fastq-view-goto-record 95)
      (should (string-match-p ":1:1:94:" (fastq-view-test--header-at-point)))
      (should (string-match-p "records 86-95 of 95" (fastq-view--header-line))))))

(ert-deftest fastq-view-test-last-page-and-count ()
  (fastq-test-with-file (f 37)
    (fastq-view-test-with (f 10)
      (fastq-view-last-page)
      (should (= fastq-view-total 37))
      (should (= fastq-record-base 27))
      (should (= (fastq-record-index-at-point) 37))
      (fastq-view-first-page)
      (should (= fastq-record-base 0)))))

(ert-deftest fastq-view-test-copy-record ()
  (fastq-view-test-with ((fastq-test-fixture "sample_R1_001.fastq") 5)
    (fastq-view-goto-record 2)
    (fastq-view-copy-record)
    (should (string-match-p "\\`@SYN01:7:FC123ABXX:1:1101:1001:2001 1:N:0:ACGTACGT\\+TTGCAACC\n[ACGTN]\\{60\\}\n\\+\n.\\{60\\}\n\\'"
                            (car kill-ring)))))

(ert-deftest fastq-view-test-gzip-file ()
  (skip-unless (or (executable-find "gzip") (fastq--zlib-p)))
  (fastq-view-test-with ((fastq-test-fixture "sample_R1_001.fastq.gz") 5)
    (fastq-view-goto-record 11)
    (should (= fastq-record-base 10))
    (should (string-match-p ":1010:2010 " (fastq-view-test--header-at-point)))
    (should (= fastq-view-total 12))))

(ert-deftest fastq-view-test-checkpoint-cache ()
  (fastq-test-with-file (f 2500)
    (let ((fastq-checkpoint-interval 1000))
      (clrhash fastq-view--checkpoints)
      (fastq-view-test-with (f 10)
        (fastq-view-goto-record 2401)
        (let ((cps (gethash (fastq-view--file-key f) fastq-view--checkpoints)))
          (should (assq 1000 cps))
          (should (assq 2000 cps))
          (should (assq 2400 cps)))
        (should (string-match-p ":1:1:2400:" (fastq-view-test--header-at-point)))
        ;; a later jump seeks from the 2000 checkpoint and lands correctly
        (fastq-view-goto-record 2051)
        (should (string-match-p ":1:1:2050:" (fastq-view-test--header-at-point)))))))

(ert-deftest fastq-view-test-reuses-buffer ()
  (let* ((f (fastq-test-fixture "sample_R1_001.fastq"))
         (a (fastq-view-noselect f))
         (b (fastq-view-noselect f 7)))
    (unwind-protect
        (progn (should (eq a b))
               (with-current-buffer b (should (= (fastq-record-index-at-point) 7))))
      (kill-buffer a))))

(ert-deftest fastq-view-test-validate-page ()
  (fastq-view-test-with ((fastq-test-fixture "invalid.fastq") 5)
    (should (string-match-p "Record 2:" (fastq-validate)))))

(provide 'fastq-view-test)
;;; fastq-view-test.el ends here
