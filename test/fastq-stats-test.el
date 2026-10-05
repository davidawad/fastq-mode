;;; fastq-stats-test.el --- Tests for fastq-stats -*- lexical-binding: t; -*-

;;; Commentary:

;; The streaming statistics are checked against a naive reference computed
;; from parsed records, on the plain and gzipped fixtures and Phred+64.

;;; Code:

(require 'fastq-test-helpers)
(require 'fastq-stats)

(defun fastq-stats-test--reference (file offset)
  "Naive statistics of FILE as an alist, qualities at OFFSET."
  (let* ((text (with-temp-buffer (insert-file-contents-literally file) (buffer-string)))
         (recs (fastq-parse-records text))
         (seqs (mapcar #'fastq-record-sequence recs))
         (quals (mapconcat #'fastq-record-quality recs ""))
         (all (apply #'concat seqs)))
    `((reads . ,(length recs))
      (bases . ,(length all))
      (gc . ,(cl-count-if (lambda (c) (memq c '(?G ?C ?g ?c))) all))
      (n . ,(cl-count-if (lambda (c) (memq c '(?N ?n))) all))
      (qsum . ,(cl-reduce #'+ (mapcar (lambda (c) (- c offset)) quals)))
      (q20 . ,(cl-count-if (lambda (c) (>= (- c offset) 20)) quals))
      (q30 . ,(cl-count-if (lambda (c) (>= (- c offset) 30)) quals)))))

(defun fastq-stats-test--check (file offset)
  "Compare `fastq-compute-stats' on FILE with the reference at OFFSET."
  (let ((st (fastq-compute-stats file))
        (ref (fastq-stats-test--reference
              (if (string-suffix-p ".gz" file) (string-remove-suffix ".gz" file) file)
              offset)))
    (should (= (fastq-stats-offset st) offset))
    (should (= (fastq-stats-reads st) (alist-get 'reads ref)))
    (should (= (fastq-stats-bases st) (alist-get 'bases ref)))
    (should (= (fastq-stats-gc st) (alist-get 'gc ref)))
    (should (= (fastq-stats-n st) (alist-get 'n ref)))
    (should (= (fastq-stats-qsum st) (alist-get 'qsum ref)))
    (should (= (fastq-stats-q20 st) (alist-get 'q20 ref)))
    (should (= (fastq-stats-q30 st) (alist-get 'q30 ref)))))

(ert-deftest fastq-stats-test-plain ()
  (fastq-stats-test--check (fastq-test-fixture "sample_R1_001.fastq") 33))

(ert-deftest fastq-stats-test-gzip ()
  (skip-unless (or (executable-find "gzip") (fastq--zlib-p)))
  (fastq-stats-test--check (fastq-test-fixture "sample_R1_001.fastq.gz") 33))

(ert-deftest fastq-stats-test-phred64 ()
  (fastq-stats-test--check (fastq-test-fixture "phred64.fastq") 64))

(ert-deftest fastq-stats-test-sample-limit ()
  (let ((st (fastq-compute-stats (fastq-test-fixture "sample_R1_001.fastq") 5)))
    (should (= (fastq-stats-reads st) 5))
    (should (= (fastq-stats-bases st) 300))))

(ert-deftest fastq-stats-test-length-distribution ()
  (with-temp-buffer
    (let ((f (make-temp-file "fastq-len-" nil ".fastq")))
      (unwind-protect
          (progn
            (let ((coding-system-for-write 'no-conversion))
              (with-temp-file f
                (insert (fastq-test-record 0 30) (fastq-test-record 1 30)
                        (fastq-test-record 2 50))))
            (let ((st (fastq-compute-stats f)))
              (should (= (fastq-stats-min-len st) 30))
              (should (= (fastq-stats-max-len st) 50))
              (should (= (gethash 30 (fastq-stats-lengths st)) 2))
              (let ((report (fastq-format-stats st)))
                (should (string-match-p "Read length     min 30, mean 36.7, max 50" report))
                (should (string-match-p "30 bp +2 reads +66.7%" report)))))
        (delete-file f)))))

(ert-deftest fastq-stats-test-report-buffer ()
  (let ((fastq-seqkit-program nil))
    (save-window-excursion
      (fastq-stats (fastq-test-fixture "sample_R1_001.fastq")))
    (let ((buf (get-buffer "*fastq stats sample_R1_001.fastq*")))
      (unwind-protect
          (with-current-buffer buf
            (should (string-match-p "Reads +12" (buffer-string)))
            (should (string-match-p "Whole file" (buffer-string)))
            (should (string-match-p "Bases >= Q30" (buffer-string))))
        (kill-buffer buf)))))

(ert-deftest fastq-stats-test-seqkit-table ()
  (should (equal (fastq--seqkit-table "file\tnum_seqs\nx.fq\t12\n")
                 "  file           x.fq\n  num_seqs       12"))
  (should (string-match-p "seqkit failed" (fastq--seqkit-table "error: boom"))))

(ert-deftest fastq-stats-test-skip-escape ()
  (should (equal (fastq--skip-escape ?^) "\\^"))
  (should (equal (fastq--at-least 30 64) "\\^-~"))
  (with-temp-buffer
    (insert "]^_~")
    (goto-char (point-min))
    (should (= (fastq--count-run-chars (point-max) (fastq--at-least 30 64)) 3))))

(provide 'fastq-stats-test)
;;; fastq-stats-test.el ends here
