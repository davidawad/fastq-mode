;;; fastq-stream-test.el --- Tests for fastq-stream -*- lexical-binding: t; -*-

;;; Commentary:

;; Skipping and taking lines from plain, gzip and zlib sources, end of file
;; handling, record counting, checkpoints and the record walker.

;;; Code:

(require 'fastq-test-helpers)
(require 'fastq-stream)

(defun fastq-stream-test--fixture-text (name)
  "Contents of fixture NAME."
  (with-temp-buffer
    (insert-file-contents-literally (fastq-test-fixture name))
    (buffer-string)))

(defun fastq-stream-test--lines (text from n)
  "Lines FROM .. FROM+N of TEXT, newline-terminated."
  (mapconcat (lambda (l) (concat l "\n"))
             (cl-subseq (split-string text "\n") from (+ from n)) ""))

(ert-deftest fastq-stream-test-gzip-magic ()
  (should (fastq-gzip-file-p (fastq-test-fixture "sample_R1_001.fastq.gz")))
  (should-not (fastq-gzip-file-p (fastq-test-fixture "sample_R1_001.fastq"))))

(ert-deftest fastq-stream-test-plain-skip-take ()
  (let* ((f (fastq-test-fixture "sample_R1_001.fastq"))
         (text (fastq-stream-test--fixture-text "sample_R1_001.fastq"))
         (res (fastq-stream-lines f 8 8)))
    (should (equal (plist-get res :text) (fastq-stream-test--lines text 8 8)))
    (should (= (plist-get res :lines) 8))
    (should (= (plist-get res :skipped) 8))
    ;; byte offset of line 8 = length of the first 8 lines
    (should (= (plist-get res :byte) (length (fastq-stream-test--lines text 0 8))))))

(ert-deftest fastq-stream-test-small-chunks ()
  "Chunks smaller than a line still give whole lines."
  (let* ((fastq-chunk-size 7)
         (f (fastq-test-fixture "sample_R1_001.fastq"))
         (text (fastq-stream-test--fixture-text "sample_R1_001.fastq")))
    (should (equal (plist-get (fastq-stream-lines f 13 9) :text)
                   (fastq-stream-test--lines text 13 9)))))

(ert-deftest fastq-stream-test-past-end ()
  (let ((res (fastq-stream-lines (fastq-test-fixture "sample_R1_001.fastq") 1000 4)))
    (should (equal (plist-get res :text) ""))
    (should (= (plist-get res :skipped) 48))
    (should (plist-get res :eof))))

(ert-deftest fastq-stream-test-no-final-newline ()
  (let ((f (fastq-test-fixture "no-final-newline.fastq")))
    (should (= (fastq-stream-count-records f) 12))
    (let ((res (fastq-stream-lines f 44 4)))
      (should (= (plist-get res :lines) 4))
      (should (string-suffix-p "\n" (plist-get res :text))))))

(ert-deftest fastq-stream-test-count ()
  (should (= (fastq-stream-count-records (fastq-test-fixture "sample_R1_001.fastq")) 12)))

(ert-deftest fastq-stream-test-zlib-source ()
  "Without a gzip executable, small .gz files go through Emacs' zlib."
  (skip-unless (fastq--zlib-p))
  (let ((fastq-gzip-program nil)
        (gz (fastq-test-fixture "sample_R1_001.fastq.gz"))
        (text (fastq-stream-test--fixture-text "sample_R1_001.fastq")))
    (should (eq (fastq-file-kind gz) 'zlib))
    (should (equal (plist-get (fastq-stream-lines gz 4 8) :text)
                   (fastq-stream-test--lines text 4 8)))
    (should (= (fastq-stream-count-records gz) 12))))

(ert-deftest fastq-stream-test-zlib-size-limit ()
  (let ((fastq-gzip-program nil) (fastq-zlib-max-size 10))
    (should-error (fastq-file-kind (fastq-test-fixture "sample_R1_001.fastq.gz"))
                  :type 'fastq-error)))

(ert-deftest fastq-stream-test-gzip-source ()
  (skip-unless (executable-find "gzip"))
  (let ((gz (fastq-test-fixture "sample_R1_001.fastq.gz"))
        (text (fastq-stream-test--fixture-text "sample_R1_001.fastq")))
    (should (eq (fastq-file-kind gz) 'gzip))
    (should (equal (plist-get (fastq-stream-lines gz 20 12) :text)
                   (fastq-stream-test--lines text 20 12)))
    (should (= (fastq-stream-count-records gz) 12))
    (should-not (cl-find-if (lambda (p) (string-prefix-p "fastq-gunzip" (process-name p)))
                            (process-list)))))

(ert-deftest fastq-stream-test-gzip-stops-early ()
  "Taking the first page of a big gzip file does not decompress all of it."
  (fastq-test-with-file (f 30000 t)
    (let ((t0 (float-time))
          (res (fastq-stream-lines f 0 8)))
      (should (= (plist-get res :lines) 8))
      (should-not (plist-get res :eof))
      (should (< (- (float-time) t0) 5)))))

(ert-deftest fastq-stream-test-checkpoints ()
  (fastq-test-with-file (f 25000)
    (let ((fastq-checkpoint-interval 10000) cps)
      (let ((res (fastq-stream-lines f (* 4 23000) 4
                                     :on-checkpoint (lambda (cp) (push cp cps)))))
        (should (string-prefix-p "@T:1:FC:1:1:23000:" (plist-get res :text))))
      (should (equal (sort (mapcar #'car cps) #'<) '(10000 20000)))
      ;; seeking from a checkpoint lands on the same record
      (let* ((cp (assq 20000 cps))
             (res (fastq-stream-lines f (* 4 3000) 4 :start-byte (cdr cp)
                                      :first-record 20000)))
        (should (string-prefix-p "@T:1:FC:1:1:23000:" (plist-get res :text)))))))

(ert-deftest fastq-stream-test-each-record ()
  (let (headers)
    (should (= (fastq-stream-each-record
                (fastq-test-fixture "no-final-newline.fastq")
                (lambda (h _s q)
                  (push (buffer-substring h (save-excursion (goto-char h) (line-end-position)))
                        headers)
                  (should (= (char-after h) ?@))
                  (should (/= (char-after q) ?+))))
               12))
    (should (string-match-p ":1011:2011 " (car headers)))))

(ert-deftest fastq-stream-test-each-record-max ()
  (let ((fastq-chunk-size 100))
    (should (= (fastq-stream-each-record (fastq-test-fixture "sample_R1_001.fastq")
                                         #'ignore 5)
               5))))

(provide 'fastq-stream-test)
;;; fastq-stream-test.el ends here
