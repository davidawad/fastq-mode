;;; screenshots.el --- Render the README screenshots as SVG/PNG -*- lexical-binding: t; -*-

;;; Commentary:

;; Run from the repository root:
;;
;;   emacs -Q --batch -l examples/screenshots.el
;;
;; Opens the synthetic demo files (examples/make-demo.py) in `fastq-view'
;; and `fastq-stats', fontifies them with the package's own faces, and
;; writes what an Emacs window would show -- header line, buffer text with
;; its faces, mode line and echo area -- as SVG in docs/screenshots/.
;; When rsvg-convert is installed, PNGs are written next to them.  No GUI
;; or X display is needed.  Colours come from each face's dark-background
;; `defface' spec; inherited standard faces use the modus-vivendi palette.

;;; Code:

(require 'cl-lib)

(defvar shot-root
  (file-name-directory
   (directory-file-name (file-name-directory (or load-file-name buffer-file-name))))
  "Repository root.")

(add-to-list 'load-path shot-root)
(require 'fastq-view)
(require 'fastq-stats)

(defconst shot-dir (expand-file-name "docs/screenshots" shot-root))

;; Show the demo files as a user's checkout, not wherever this runs.
(setq directory-abbrev-alist
      (list (cons (concat "\\`" (regexp-quote shot-root)) "~/src/fastq-mode/")))

(defconst shot-palette
  '((bg . "#0d0e1c") (fg . "#ffffff") (bar . "#2d2f42") (bar-fg . "#e0e0e0")
    (font-lock-function-name-face . "#feacd0")
    (font-lock-comment-face . "#989898")
    (shadow . "#989898"))
  "Colours for the frame and for standard faces the package inherits.")

(defun shot--face-fg (face)
  "Foreground colour of FACE for a dark background."
  (let* ((spec (get face 'face-defface-spec))
         (entry (or (cl-find-if (lambda (e) (equal (car e) '((background dark)))) spec)
                    (assq t spec)))
         (atts (cdr entry)))
    (or (plist-get atts :foreground)
        (alist-get (plist-get atts :inherit) shot-palette)
        (alist-get face shot-palette)
        (alist-get 'fg shot-palette))))

(defun shot--style (faces)
  "SVG attributes for the list of FACES."
  (let* ((faces (if (listp faces) faces (list faces)))
         (main (cl-find-if (lambda (f) (not (eq f 'fastq-low-quality-base-face))) faces))
         (low (memq 'fastq-low-quality-base-face faces)))
    (concat (format " fill=\"%s\"" (if main (shot--face-fg main) (alist-get 'fg shot-palette)))
            (if (memq main '(fastq-base-n-face fastq-header-face)) " font-weight=\"bold\"" "")
            (if low " text-decoration=\"underline\" fill-opacity=\"0.55\"" ""))))

(defun shot--escape (s)
  "S escaped for XML."
  (replace-regexp-in-string
   "[<>&]" (lambda (m) (pcase m ("<" "&lt;") (">" "&gt;") (_ "&amp;"))) s t t))

(defun shot--line-svg (beg end)
  "SVG tspans for the buffer text BEG..END, split where faces change."
  (let ((pos beg) (out ""))
    (while (< pos end)
      (let* ((next (min end (next-single-property-change pos 'face nil end)))
             (text (buffer-substring-no-properties pos next)))
        (setq out (concat out (format "<tspan%s>%s</tspan>"
                                      (shot--style (get-text-property pos 'face))
                                      (shot--escape text))))
        (setq pos next)))
    out))

(defun shot--bar (y w text cw)
  "A header/mode-line bar at Y, W wide, with TEXT; CW is the char width."
  (ignore cw)
  (format "<rect x=\"0\" y=\"%d\" width=\"%d\" height=\"22\" fill=\"%s\"/><text x=\"12\" y=\"%d\" fill=\"%s\" xml:space=\"preserve\">%s</text>\n"
          y w (alist-get 'bar shot-palette) (+ y 16) (alist-get 'bar-fg shot-palette)
          (shot--escape text)))

(defun shot-render (name &rest args)
  "Write docs/screenshots/NAME.svg (and .png) from the current buffer.
ARGS: :lines N, :cols N, :header STRING, :mode STRING, :echo STRING."
  (let* ((lines (or (plist-get args :lines) 30)) (cols (or (plist-get args :cols) 112))
         (cw 8.43) (lh 19) (w (round (+ 24 (* cw cols))))
         (top (if (plist-get args :header) 30 8))
         (h (+ top (* lh lines) 60))
         (svg (list (format "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" font-family=\"JetBrains Mono, Menlo, DejaVu Sans Mono, monospace\" font-size=\"14\">\n<rect width=\"100%%\" height=\"100%%\" fill=\"%s\"/>\n"
                            w h (alist-get 'bg shot-palette)))))
    (when (plist-get args :header)
      (push (shot--bar 0 w (plist-get args :header) cw) svg))
    (save-excursion
      (goto-char (point-min))
      (dotimes (i lines)
        (unless (eobp)
          (let ((end (min (line-end-position) (+ (point) cols))))
            (push (format "<text x=\"12\" y=\"%d\" xml:space=\"preserve\">%s</text>\n"
                          (+ top 15 (* i lh)) (shot--line-svg (point) end))
                  svg))
          (forward-line 1))))
    (push (shot--bar (+ top (* lh lines) 6) w (or (plist-get args :mode) "") cw) svg)
    (when (plist-get args :echo)
      (push (format "<text x=\"12\" y=\"%d\" fill=\"%s\" xml:space=\"preserve\">%s</text>\n"
                    (+ top (* lh lines) 48) (alist-get 'fg shot-palette)
                    (shot--escape (plist-get args :echo)))
            svg))
    (push "</svg>\n" svg)
    (make-directory shot-dir t)
    (let ((file (expand-file-name (concat name ".svg") shot-dir)))
      (with-temp-file file (insert (apply #'concat (nreverse svg))))
      (let ((rsvg (executable-find "rsvg-convert")))
        (when rsvg
          (call-process rsvg nil nil nil "-z" "1.5" "-o"
                        (expand-file-name (concat name ".png") shot-dir) file)))
      (message "wrote %s" file))))

;;;; The shots

(let* ((demo (expand-file-name "examples/demo_R1_001.fastq.gz" shot-root))
       (fastq-view-page-records 200)
       (buf (fastq-view-noselect demo 41)))
  (with-current-buffer buf
    (font-lock-ensure)
    ;; `format-mode-line' needs a window, which batch Emacs lacks
    (let ((header (fastq-view--header-line)))
      (goto-char (point-min))
      ;; put point on a weak base so the echo area shows eldoc for it
      (forward-line (+ 4 1))
      (let* ((qual (save-excursion (forward-line 2)
                                   (buffer-substring (point) (line-end-position))))
             (col (or (string-match "[#-]" qual) 10)))
        (forward-char col)
        (let ((echo (fastq-eldoc-function)))
          (save-restriction
            (goto-char (point-min))
            (forward-line 0)
            (narrow-to-region (point) (point-max))
            (shot-render "viewer" :lines 28 :header header
                         :mode " *fastq demo_R1_001.fastq.gz*   (FASTQ-view)"
                         :echo echo)))))))

(let ((fastq-seqkit-program nil))
  (save-window-excursion
    (fastq-stats (expand-file-name "examples/demo_R1_001.fastq.gz" shot-root)))
  (with-current-buffer "*fastq stats demo_R1_001.fastq.gz*"
    (goto-char (point-min))
    ;; drop the timing line, which changes on every run
    (let ((inhibit-read-only t))
      (when (re-search-forward "^(.* s)\n" nil t) (replace-match "")))
    (shot-render "stats" :lines 24 :cols 72
                 :mode " *fastq stats demo_R1_001.fastq.gz*   (Special)")))

(provide 'screenshots)
;;; screenshots.el ends here
