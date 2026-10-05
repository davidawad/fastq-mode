;;; fastq-mode-test.el --- Tests for fastq-mode, highlighting and eldoc -*- lexical-binding: t; -*-

;;; Commentary:

;; The major mode on small files: record-aware faces, eldoc, navigation,
;; validation, header help, mate files and `fastq-open' dispatch.

;;; Code:

(require 'fastq-test-helpers)
(require 'fastq-mode)
(require 'fastq-view)

(defmacro fastq-mode-test-with-fixture (name &rest body)
  "Visit fixture NAME in `fastq-mode', fontified, and run BODY."
  (declare (indent 1))
  `(let ((buf (find-file-noselect (fastq-test-fixture ,name))))
     (unwind-protect
         (with-current-buffer buf
           (fastq-mode)
           (font-lock-ensure)
           (goto-char (point-min))
           ,@body)
       (kill-buffer buf))))

(ert-deftest fastq-mode-test-auto-mode ()
  (should (eq (cdr (assoc "\\.f\\(?:ast\\)?q\\'" auto-mode-alist)) 'fastq-mode))
  (let ((buf (find-file-noselect (fastq-test-fixture "sample_R1_001.fastq"))))
    (unwind-protect
        (should (eq (buffer-local-value 'major-mode buf) 'fastq-mode))
      (kill-buffer buf))))

(ert-deftest fastq-mode-test-faces-by-role ()
  (fastq-mode-test-with-fixture "sample_R1_001.fastq"
    (should (memq 'fastq-header-face (fastq-test-faces-at (point))))
    (should (memq 'fastq-header-comment-face (fastq-test-faces-at (1- (line-end-position)))))
    (forward-line 1)
    (let ((base (char-after)))
      (should (memq (fastq-base-face base) (fastq-test-faces-at (point)))))
    (forward-line 1)
    (should (memq 'fastq-separator-face (fastq-test-faces-at (point))))
    (forward-line 1)
    (should (memq (fastq-quality-face (char-after) 33) (fastq-test-faces-at (point))))))

(ert-deftest fastq-mode-test-quality-line-starting-with-at ()
  "A quality line that begins with @ is still coloured as qualities."
  (with-temp-buffer
    (insert "@r1\nACGT\n+\n@@II\n@r2\nACGT\n+\nIIII\n")
    (fastq-mode)
    (font-lock-ensure)
    (goto-char (point-min))
    (forward-line 3)
    (should (memq 'fastq-quality-good-face (fastq-test-faces-at (point))))
    (should-not (memq 'fastq-header-face (fastq-test-faces-at (point))))))

(ert-deftest fastq-mode-test-low-quality-bases-marked ()
  (with-temp-buffer
    (insert "@r1\nACGT\n+\n#III\n")
    (fastq-mode)
    (font-lock-ensure)
    (goto-char (point-min))
    (forward-line 1)
    (should (memq 'fastq-low-quality-base-face (fastq-test-faces-at (point))))
    (should-not (memq 'fastq-low-quality-base-face (fastq-test-faces-at (1+ (point)))))
    (let ((fastq-shade-low-quality-bases nil))
      (font-lock-flush) (font-lock-ensure)
      (should-not (memq 'fastq-low-quality-base-face (fastq-test-faces-at (point)))))))

(ert-deftest fastq-mode-test-phred64-detected ()
  (fastq-mode-test-with-fixture "phred64.fastq"
    (should (= fastq-buffer-offset 64))))

(ert-deftest fastq-mode-test-eldoc ()
  (with-temp-buffer
    (insert "@A:1:FC:2:3:4:5 1:N:0:ACGT\nACGT\n+\nI#5I\n")
    (fastq-mode)
    (goto-char (point-min))
    (should (string-match-p "Read 1: Illumina CASAVA 1.8\\+.*lane 2" (fastq-eldoc-function)))
    (forward-line 1)
    (forward-char 1)
    (should (string-match-p "Read 1, base 2/4: C  Q2 (error 63.1%, low)" (fastq-eldoc-function)))
    (forward-line 2)
    (forward-char 2)
    (should (string-match-p "base 3/4: G  Q20 (error 1%, fair)" (fastq-eldoc-function)))
    (end-of-line)
    (should (string-match-p "Read 1: 4 bases, mean Q" (fastq-eldoc-function)))))

(ert-deftest fastq-mode-test-navigation ()
  (fastq-mode-test-with-fixture "sample_R1_001.fastq"
    (forward-line 2)
    (fastq-next-record)
    (should (= (line-number-at-pos) 5))
    (should (= (fastq-record-index-at-point) 2))
    (fastq-next-record 3)
    (should (= (fastq-record-index-at-point) 5))
    (fastq-previous-record)
    (should (= (fastq-record-index-at-point) 4))
    (fastq-goto-record 12)
    (should (looking-at "@SYN01:7:FC123ABXX:1:1101:1011:2011 "))))

(ert-deftest fastq-mode-test-record-at-point ()
  (fastq-mode-test-with-fixture "sample_R1_001.fastq"
    (forward-line 6)
    (let ((rec (fastq-record-at-point)))
      (should (string-match-p ":1001:2001 " (fastq-record-header rec)))
      (should (= (length (fastq-record-sequence rec)) 60))
      (should (= (length (fastq-record-quality rec)) 60)))))

(ert-deftest fastq-mode-test-validate ()
  (fastq-mode-test-with-fixture "sample_R1_001.fastq"
    (should (string-match-p "12 records, all valid" (fastq-validate))))
  (fastq-mode-test-with-fixture "invalid.fastq"
    (should (string-match-p "Record 2: sequence has 4 bases but quality has 3"
                            (fastq-validate)))
    (should (= (fastq-record-index-at-point) 2))))

(ert-deftest fastq-mode-test-describe-header ()
  (fastq-mode-test-with-fixture "sample_R1_001.fastq"
    (fastq-describe-header)
    (unwind-protect
        (with-current-buffer "*fastq header*"
          (should (string-match-p "flowcell +FC123ABXX" (buffer-string)))
          (should (string-match-p "index +ACGTACGT\\+TTGCAACC" (buffer-string))))
      (kill-buffer "*fastq header*"))))

(ert-deftest fastq-mode-test-mate ()
  (fastq-mode-test-with-fixture "sample_R1_001.fastq"
    (fastq-goto-record 3)
    (save-window-excursion
      (fastq-mate)
      (unwind-protect
          (progn
            (should (string-suffix-p "sample_R2_001.fastq" buffer-file-name))
            (should (= (fastq-record-index-at-point) 3))
            (should (looking-at "@SYN01:7:FC123ABXX:1:1101:1002:2002 2:")))
        (kill-buffer)))))

(ert-deftest fastq-mode-test-open-dispatch ()
  (let ((small (fastq-open-noselect (fastq-test-fixture "sample_R1_001.fastq") 2)))
    (unwind-protect
        (with-current-buffer small
          (should (eq major-mode 'fastq-mode))
          (should (= (fastq-record-index-at-point) 2)))
      (kill-buffer small)))
  (skip-unless (or (executable-find "gzip") (fastq--zlib-p)))
  (let ((gz (fastq-open-noselect (fastq-test-fixture "sample_R1_001.fastq.gz"))))
    (unwind-protect
        (with-current-buffer gz
          (should (eq major-mode 'fastq-view-mode))
          (should-not buffer-file-name))
      (kill-buffer gz))))

(ert-deftest fastq-mode-test-auto-view-mode ()
  (skip-unless (or (executable-find "gzip") (fastq--zlib-p)))
  (unwind-protect
      (progn
        (fastq-auto-view-mode 1)
        (let ((buf (find-file-noselect (fastq-test-fixture "sample_R1_001.fastq.gz"))))
          (unwind-protect
              (with-current-buffer buf
                (should (eq major-mode 'fastq-view-mode))
                (should (looking-at "@SYN01")))
            (kill-buffer buf)))
        (let* ((fastq-view-threshold 100)
               (buf (find-file-noselect (fastq-test-fixture "sample_R2_001.fastq"))))
          (unwind-protect
              (should (eq (buffer-local-value 'major-mode buf) 'fastq-view-mode))
            (kill-buffer buf))))
    (fastq-auto-view-mode -1))
  (should-not (advice-member-p #'fastq--find-file-advice 'find-file-noselect)))

(provide 'fastq-mode-test)
;;; fastq-mode-test.el ends here
