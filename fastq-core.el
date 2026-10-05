;;; fastq-core.el --- Records, Phred qualities and read headers for fastq-mode -*- lexical-binding: t; -*-

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

;; Pure functions shared by the rest of fastq-mode: the four-line record
;; layout, Phred quality decoding (offset 33 or 64), quality bins, and the
;; Illumina read-header decoder (CASAVA 1.8+ and the older pre-1.8 form).
;; Nothing here touches files or processes.

;;; Code:

(require 'cl-lib)
(require 'subr-x)

(defgroup fastq nil
  "View and inspect FASTQ sequencing reads."
  :group 'data
  :prefix "fastq-")

(define-error 'fastq-error "FASTQ error")

(defun fastq--error (fmt &rest args)
  "Signal `fastq-error' with a message built from FMT and ARGS."
  (signal 'fastq-error (list (apply #'format-message fmt args))))

;;;; Record layout

;; FASTQ as written by every current instrument and converter: four lines
;; per record (@header, sequence, +[header], qualities), no line wrapping.
;; Wrapped (multi-line) FASTQ is rejected by `fastq-validate-record'.

(defconst fastq-lines-per-record 4
  "Lines per FASTQ record.")

(defun fastq-line-role (line)
  "Return the role of 1-based LINE, counting from a record boundary.
One of `header', `sequence', `separator' or `quality'."
  (aref [header sequence separator quality]
        (mod (1- line) fastq-lines-per-record)))

(cl-defstruct (fastq-record (:constructor fastq-record-create)
                            (:copier nil))
  "One FASTQ read."
  header sequence quality)

(defun fastq-parse-records (text)
  "Parse TEXT (whole four-line records) into a list of `fastq-record'.
A trailing partial record is ignored."
  (let ((lines (split-string text "\n"))
        out)
    (while (nthcdr 3 lines)
      (let ((h (pop lines)) (s (pop lines)) (_ (pop lines)) (q (pop lines)))
        (push (fastq-record-create :header h :sequence s :quality q) out)))
    (nreverse out)))

(defun fastq-validate-record (header sequence separator quality)
  "Return nil if HEADER SEQUENCE SEPARATOR QUALITY form a valid record.
Otherwise return a string describing the first problem."
  (cond
   ((not (string-prefix-p "@" header)) "header line does not start with @")
   ((not (string-prefix-p "+" separator)) "third line does not start with +")
   ((/= (length sequence) (length quality))
    (format "sequence has %d bases but quality has %d characters"
            (length sequence) (length quality)))
   ((string-match-p "[^ACGTNacgtnRYSWKMBDHVryswkmbdhv.-]" sequence)
    "sequence contains characters that are not IUPAC bases")))

;;;; Phred qualities

(defcustom fastq-phred-offset 'auto
  "ASCII offset of quality characters: 33, 64 or `auto' (detect per file).
Every instrument since 2011 writes Phred+33 (Sanger / Illumina 1.8+)."
  :type '(choice (const :tag "Detect" auto) (const 33) (const 64)))

(defun fastq-detect-offset (qualities)
  "Guess the Phred offset from QUALITIES, a list of quality strings.
Return 33 or 64.  Any character below `;' (59) can only be Phred+33;
an alphabet that stays at or above `@' (64) and reaches past `J' (74)
is the old Illumina 1.3-1.7 Phred+64.  Ambiguous input is 33."
  (let ((lo 255) (hi 0))
    (dolist (q qualities)
      (dotimes (i (length q))
        (let ((c (aref q i)))
          (setq lo (min lo c) hi (max hi c)))))
    (if (and (>= lo 64) (> hi 74)) 64 33)))

(defun fastq-offset (qualities)
  "The Phred offset to use for QUALITIES, honouring `fastq-phred-offset'."
  (if (eq fastq-phred-offset 'auto)
      (fastq-detect-offset qualities)
    fastq-phred-offset))

(defun fastq-phred (char offset)
  "Phred score encoded by quality CHAR with OFFSET."
  (- char offset))

(defun fastq-error-probability (q)
  "Probability that a base called with Phred score Q is wrong."
  (expt 10.0 (/ (- q) 10.0)))

(defcustom fastq-quality-bins '(10 20 30)
  "Phred thresholds separating the low, poor, fair and good quality bins."
  :type '(list integer integer integer))

(defun fastq-quality-bin (q)
  "Bin of Phred score Q: `low', `poor', `fair' or `good'."
  (pcase-let ((`(,a ,b ,c) fastq-quality-bins))
    (cond ((< q a) 'low) ((< q b) 'poor) ((< q c) 'fair) (t 'good))))

(defun fastq-mean-quality (quality offset)
  "Mean Phred score of the QUALITY string with OFFSET (0.0 when empty)."
  (if (string-empty-p quality)
      0.0
    (/ (float (- (cl-reduce #'+ quality) (* offset (length quality))))
       (length quality))))

;;;; Read headers

(defun fastq--casava18 (id comment)
  "Decode an Illumina CASAVA 1.8+ header from ID and COMMENT, or nil.
ID is @instrument:run:flowcell:lane:tile:x:y[:umi], COMMENT is
read:filtered:control:index."
  (let ((f (split-string id ":")))
    (when (and (memq (length f) '(7 8))
               (cl-every (lambda (s) (string-match-p "\\`[0-9]+\\'" s))
                         (cl-subseq f 3 7)))
      (let ((c (and comment (split-string comment ":"))))
        (append
         `((format . "Illumina CASAVA 1.8+")
           (instrument . ,(nth 0 f)) (run . ,(nth 1 f)) (flowcell . ,(nth 2 f))
           (lane . ,(nth 3 f)) (tile . ,(nth 4 f)) (x . ,(nth 5 f)) (y . ,(nth 6 f)))
         (and (nth 7 f) `((umi . ,(nth 7 f))))
         (and (= (length c) 4)
              `((read . ,(nth 0 c))
                (filtered . ,(if (equal (nth 1 c) "Y") "yes (failed filter)" "no"))
                (control . ,(nth 2 c))
                (index . ,(nth 3 c)))))))))

(defun fastq--casava-old (id)
  "Decode a pre-1.8 Illumina ID (instrument:lane:tile:x:y#index/read), or nil."
  (when (string-match
         "\\`\\([^:]+\\):\\([0-9]+\\):\\([0-9]+\\):\\([0-9]+\\):\\([0-9]+\\)\\(?:#\\([^/]*\\)\\)?\\(?:/\\([12]\\)\\)?\\'"
         id)
    (let ((m (lambda (n) (match-string n id))))
      (cl-remove-if-not
       #'cdr
       `((format . "Illumina pre-1.8")
         (instrument . ,(funcall m 1)) (lane . ,(funcall m 2)) (tile . ,(funcall m 3))
         (x . ,(funcall m 4)) (y . ,(funcall m 5))
         (index . ,(funcall m 6)) (read . ,(funcall m 7)))))))

(defun fastq-split-header (header)
  "Split HEADER (with or without @) into (ID . COMMENT); COMMENT may be nil."
  (let* ((h (string-remove-prefix "@" header))
         (sp (string-search " " h)))
    (if sp (cons (substring h 0 sp) (substring h (1+ sp))) (cons h nil))))

(defun fastq-decode-header (header)
  "Decode read HEADER into an alist of named fields.
Always contains `id'; Illumina headers also carry instrument, run,
flowcell, lane, tile, x, y and (CASAVA 1.8+) read, filtered, control
and index.  Unknown formats return just `id' and `comment'."
  (pcase-let ((`(,id . ,comment) (fastq-split-header header)))
    (append `((id . ,id))
            (or (fastq--casava18 id comment)
                (fastq--casava-old id)
                (and comment `((comment . ,comment)))))))

(defun fastq-read-name (header)
  "The read name shared by both mates of HEADER (no @, comment or /1 /2)."
  (replace-regexp-in-string "/[12]\\'" "" (car (fastq-split-header header))))

(defun fastq-header-summary (header)
  "One-line human summary of HEADER."
  (let ((d (fastq-decode-header header)))
    (if (alist-get 'lane d)
        (string-join
         (delq nil
               (list (format "%s" (alist-get 'format d))
                     (format "instrument %s" (alist-get 'instrument d))
                     (and (alist-get 'flowcell d)
                          (format "flowcell %s" (alist-get 'flowcell d)))
                     (format "lane %s tile %s (x %s, y %s)" (alist-get 'lane d)
                             (alist-get 'tile d) (alist-get 'x d) (alist-get 'y d))
                     (and (alist-get 'read d) (format "read %s" (alist-get 'read d)))
                     (and (alist-get 'index d) (format "index %s" (alist-get 'index d)))
                     (and (alist-get 'filtered d)
                          (format "filtered %s" (alist-get 'filtered d)))))
         ", ")
      (format "read %s" (alist-get 'id d)))))

;;;; Mates

(defconst fastq--mate-patterns
  '(("_R1_" . "_R2_") ("_R1." . "_R2.") ("_1." . "_2.") (".R1." . ".R2.")
    ("_r1_" . "_r2_") ("_read1" . "_read2"))
  "Pairs of file-name fragments that distinguish read 1 from read 2.")

(defun fastq-mate-file-name (file)
  "Name of the paired-end mate of FILE (R1 <-> R2), or nil if none is implied.
Only the last matching fragment in the base name is swapped."
  (let ((dir (file-name-directory file))
        (base (file-name-nondirectory file)))
    (cl-some
     (lambda (pair)
       (cl-some
        (lambda (from-to)
          (let ((i (string-search (car from-to) base)))
            (when i
              ;; prefer the last occurrence (sample names may contain _1.)
              (let ((last i))
                (while (setq i (string-search (car from-to) base (1+ i)))
                  (setq last i))
                (concat dir (substring base 0 last) (cdr from-to)
                        (substring base (+ last (length (car from-to)))))))))
        (list pair (cons (cdr pair) (car pair)))))
     fastq--mate-patterns)))

(provide 'fastq-core)
;;; fastq-core.el ends here
