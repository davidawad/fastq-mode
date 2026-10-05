;;; fastq-font.el --- Highlighting for FASTQ buffers -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Keywords: data, tools, faces
;; URL: https://github.com/davidawad/fastq-mode

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; FASTQ cannot be highlighted with regexps alone: a quality line may start
;; with `@' or `+', so a line's role comes from its position in the
;; four-line record.  `fastq-fontify-region' (run by font-lock through a
;; function matcher) works out each line's role from its line number, colours bases by
;; nucleotide, quality characters by Phred bin, and, when
;; `fastq-shade-low-quality-bases' is on, dims sequence bases whose call
;; is below `fastq-low-quality-threshold' -- the quality line read two
;; lines further down, so the alignment is exact.

;;; Code:

(require 'fastq-core)

(defgroup fastq-faces nil
  "Faces used by fastq-mode."
  :group 'fastq
  :group 'faces)

(defface fastq-header-face '((t :inherit font-lock-function-name-face))
  "Read header line (@...)."
  :group 'fastq-faces)
(defface fastq-header-comment-face '((t :inherit font-lock-comment-face))
  "Comment part of a read header (after the first space)."
  :group 'fastq-faces)
(defface fastq-separator-face '((t :inherit shadow))
  "The + separator line."
  :group 'fastq-faces)

(defface fastq-base-a-face '((((background dark)) :foreground "#5fd75f")
                             (t :foreground "#008700"))
  "Adenine."
  :group 'fastq-faces)
(defface fastq-base-c-face '((((background dark)) :foreground "#5fafff")
                             (t :foreground "#0000d7"))
  "Cytosine."
  :group 'fastq-faces)
(defface fastq-base-g-face '((((background dark)) :foreground "#ffd75f")
                             (t :foreground "#af8700"))
  "Guanine."
  :group 'fastq-faces)
(defface fastq-base-t-face '((((background dark)) :foreground "#ff5f5f")
                             (t :foreground "#d70000"))
  "Thymine."
  :group 'fastq-faces)
(defface fastq-base-n-face '((t :inherit shadow :weight bold))
  "N (no call) and other IUPAC codes."
  :group 'fastq-faces)
(defface fastq-low-quality-base-face '((t :underline (:style wave :color "gray50")))
  "Added to bases whose quality is below `fastq-low-quality-threshold'."
  :group 'fastq-faces)

(defface fastq-quality-low-face '((((background dark)) :foreground "#ff5f5f")
                                  (t :foreground "#d70000"))
  "Quality characters in the lowest bin (see `fastq-quality-bins')."
  :group 'fastq-faces)
(defface fastq-quality-poor-face '((((background dark)) :foreground "#ffaf5f")
                                   (t :foreground "#d75f00"))
  "Quality characters in the second bin."
  :group 'fastq-faces)
(defface fastq-quality-fair-face '((((background dark)) :foreground "#d7d75f")
                                   (t :foreground "#878700"))
  "Quality characters in the third bin."
  :group 'fastq-faces)
(defface fastq-quality-good-face '((((background dark)) :foreground "#5faf5f")
                                   (t :foreground "#005f00"))
  "Quality characters at or above the top threshold."
  :group 'fastq-faces)

(defcustom fastq-shade-low-quality-bases t
  "Non-nil means mark sequence bases whose quality is low."
  :type 'boolean
  :group 'fastq)

(defcustom fastq-low-quality-threshold 20
  "Phred score below which a base is marked as low quality."
  :type 'integer
  :group 'fastq)

(defvar-local fastq-buffer-offset 33
  "Phred offset of the qualities in this buffer (33 or 64).")

(defun fastq-base-face (char)
  "Face for base CHAR."
  (pcase (upcase char)
    (?A 'fastq-base-a-face) (?C 'fastq-base-c-face)
    (?G 'fastq-base-g-face) (?T 'fastq-base-t-face)
    (_ 'fastq-base-n-face)))

(defun fastq-quality-face (char offset)
  "Face for quality CHAR with OFFSET."
  (pcase (fastq-quality-bin (fastq-phred char offset))
    ('low 'fastq-quality-low-face) ('poor 'fastq-quality-poor-face)
    ('fair 'fastq-quality-fair-face) (_ 'fastq-quality-good-face)))

(defun fastq--put-runs (beg end face-fn)
  "Apply faces from FACE-FN (index -> face) to BEG..END, run by run."
  (let ((i beg))
    (while (< i end)
      (let* ((face (funcall face-fn (- i beg))) (j (1+ i)))
        (while (and (< j end) (equal (funcall face-fn (- j beg)) face))
          (setq j (1+ j)))
        (put-text-property i j 'face face)
        (setq i j)))))

(defun fastq--quality-below (beg)
  "Quality string two lines below the line at BEG, or nil."
  (save-excursion
    (goto-char beg)
    (and (zerop (forward-line 2))
         (buffer-substring-no-properties (point) (line-end-position)))))

(defun fastq--fontify-sequence (beg end)
  "Highlight the sequence line BEG..END."
  (let* ((seq (buffer-substring-no-properties beg end))
         (qual (and fastq-shade-low-quality-bases (fastq--quality-below beg)))
         (off fastq-buffer-offset))
    (fastq--put-runs
     beg end
     (lambda (i)
       (let ((face (fastq-base-face (aref seq i))))
         (if (and qual (< i (length qual))
                  (< (fastq-phred (aref qual i) off) fastq-low-quality-threshold))
             (list 'fastq-low-quality-base-face face)
           face))))))

(defun fastq--fontify-quality (beg end)
  "Highlight the quality line BEG..END."
  (let ((qual (buffer-substring-no-properties beg end)) (off fastq-buffer-offset))
    (fastq--put-runs beg end (lambda (i) (fastq-quality-face (aref qual i) off)))))

(defun fastq--fontify-header (beg end)
  "Highlight the header line BEG..END."
  (let ((sp (save-excursion (goto-char beg) (search-forward " " end t))))
    (put-text-property beg (or sp end) 'face 'fastq-header-face)
    (when sp (put-text-property sp end 'face 'fastq-header-comment-face))))

(defun fastq-fontify-region (beg end)
  "Highlight whole lines from BEG to END by their role in the record."
  (save-excursion
    (with-silent-modifications
      (goto-char beg)
      (forward-line 0)
      (let ((line (line-number-at-pos)))
        (while (and (< (point) end) (not (eobp)))
          (let ((b (point)) (e (line-end-position)))
            (remove-text-properties b e '(face nil))
            (pcase (fastq-line-role line)
              ('header (fastq--fontify-header b e))
              ('sequence (fastq--fontify-sequence b e))
              ('separator (put-text-property b e 'face 'fastq-separator-face))
              ('quality (fastq--fontify-quality b e)))
            (forward-line 1)
            (setq line (1+ line))))))))

(defun fastq--font-lock-matcher (limit)
  "Font-lock matcher: highlight from point to LIMIT, never report a match."
  (fastq-fontify-region (point) limit)
  (goto-char limit)
  nil)

(defconst fastq-font-lock-keywords '((fastq--font-lock-matcher))
  "Font-lock keywords for FASTQ: a single function doing record-aware work.")

(defun fastq-font-setup ()
  "Highlight the current buffer with record-aware font-lock."
  (setq-local font-lock-defaults '(fastq-font-lock-keywords t)))

(provide 'fastq-font)
;;; fastq-font.el ends here
