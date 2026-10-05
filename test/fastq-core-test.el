;;; fastq-core-test.el --- Tests for fastq-core -*- lexical-binding: t; -*-

;;; Commentary:

;; Record parsing, Phred offsets and bins, header decoding, mate names.

;;; Code:

(require 'fastq-test-helpers)
(require 'fastq-core)

(ert-deftest fastq-core-test-line-roles ()
  (should (equal (mapcar #'fastq-line-role '(1 2 3 4 5 8))
                 '(header sequence separator quality header quality))))

(ert-deftest fastq-core-test-parse-records ()
  (let ((recs (fastq-parse-records "@a\nAC\n+\nII\n@b\nGT\n+\n##\n@partial\nA\n")))
    (should (= (length recs) 2))
    (should (equal (fastq-record-header (cadr recs)) "@b"))
    (should (equal (fastq-record-quality (car recs)) "II"))))

(ert-deftest fastq-core-test-validate ()
  (should-not (fastq-validate-record "@r" "ACGTN" "+" "IIIII"))
  (should (string-match-p "does not start with @" (fastq-validate-record "r" "A" "+" "I")))
  (should (string-match-p "third line" (fastq-validate-record "@r" "A" "-" "I")))
  (should (string-match-p "4 bases but quality has 3"
                          (fastq-validate-record "@r" "ACGT" "+" "III")))
  (should (string-match-p "IUPAC" (fastq-validate-record "@r" "ACXT" "+" "IIII"))))

(ert-deftest fastq-core-test-offset-detection ()
  (should (= (fastq-detect-offset '("II9-#")) 33))
  (should (= (fastq-detect-offset '("hhhJJ@")) 64))
  ;; ambiguous (all within 64..74) defaults to 33
  (should (= (fastq-detect-offset '("@ABCDEFGHIJ")) 33))
  (let ((fastq-phred-offset 64))
    (should (= (fastq-offset '("!!!")) 64))))

(ert-deftest fastq-core-test-phred-math ()
  (should (= (fastq-phred ?I 33) 40))
  (should (< (abs (- (fastq-error-probability 30) 0.001)) 1e-12))
  (should (eq (fastq-quality-bin 5) 'low))
  (should (eq (fastq-quality-bin 15) 'poor))
  (should (eq (fastq-quality-bin 25) 'fair))
  (should (eq (fastq-quality-bin 30) 'good))
  (should (= (fastq-mean-quality "I#" 33) 21.0))
  (should (= (fastq-mean-quality "" 33) 0.0)))

(ert-deftest fastq-core-test-casava18-header ()
  (let ((d (fastq-decode-header
            "@LH00999:26:22DEMOLT4:1:1101:12345:1000 1:N:0:ACGTGACTAG+TTAGCCATGA")))
    (should (equal (alist-get 'instrument d) "LH00999"))
    (should (equal (alist-get 'run d) "26"))
    (should (equal (alist-get 'flowcell d) "22DEMOLT4"))
    (should (equal (alist-get 'lane d) "1"))
    (should (equal (alist-get 'tile d) "1101"))
    (should (equal (alist-get 'x d) "12345"))
    (should (equal (alist-get 'y d) "1000"))
    (should (equal (alist-get 'read d) "1"))
    (should (equal (alist-get 'filtered d) "no"))
    (should (equal (alist-get 'index d) "ACGTGACTAG+TTAGCCATGA"))))

(ert-deftest fastq-core-test-casava18-umi-and-filtered ()
  (let ((d (fastq-decode-header "@M1:5:FC:2:11:3:4:ACGTTT 2:Y:0:7")))
    (should (equal (alist-get 'umi d) "ACGTTT"))
    (should (equal (alist-get 'filtered d) "yes (failed filter)"))))

(ert-deftest fastq-core-test-old-illumina-header ()
  (let ((d (fastq-decode-header "@HWUSI-EAS100R:6:73:941:1973#0/1")))
    (should (equal (alist-get 'format d) "Illumina pre-1.8"))
    (should (equal (alist-get 'lane d) "6"))
    (should (equal (alist-get 'index d) "0"))
    (should (equal (alist-get 'read d) "1"))))

(ert-deftest fastq-core-test-unknown-header ()
  (let ((d (fastq-decode-header "@SRR001666.1 071112_SLXA length=36")))
    (should (equal (alist-get 'id d) "SRR001666.1"))
    (should (equal (alist-get 'comment d) "071112_SLXA length=36"))
    (should (string-match-p "read SRR001666.1"
                            (fastq-header-summary "@SRR001666.1 x")))))

(ert-deftest fastq-core-test-header-summary ()
  (should (string-match-p
           "lane 1 tile 1101 (x 12345, y 1000), read 1"
           (fastq-header-summary "@LH00999:26:22DEMOLT4:1:1101:12345:1000 1:N:0:AC"))))

(ert-deftest fastq-core-test-read-name ()
  (should (equal (fastq-read-name "@A:1:B:1:1:2:3 1:N:0:X") "A:1:B:1:1:2:3"))
  (should (equal (fastq-read-name "@HW:6:73:941:1973#0/2") "HW:6:73:941:1973#0")))

(ert-deftest fastq-core-test-mate-file-name ()
  (let ((dir (file-name-as-directory (expand-file-name "x" temporary-file-directory))))
    (should (equal (fastq-mate-file-name (concat dir "S_L001_R1_001.fastq.gz"))
                   (concat dir "S_L001_R2_001.fastq.gz")))
    (should (equal (fastq-mate-file-name (concat dir "S_L001_R2_001.fastq.gz"))
                   (concat dir "S_L001_R1_001.fastq.gz")))
    (should (equal (fastq-mate-file-name (concat dir "run_1.fq")) (concat dir "run_2.fq")))
    ;; the last occurrence is swapped, not one inside the sample name
    (should (equal (fastq-mate-file-name (concat dir "a_1.b_1.fq")) (concat dir "a_1.b_2.fq")))
    (should-not (fastq-mate-file-name (concat dir "single.fastq")))))

(provide 'fastq-core-test)
;;; fastq-core-test.el ends here
