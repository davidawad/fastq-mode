;;; run-tests.el --- Batch entry point for tests, byte-compile and checkdoc -*- lexical-binding: t; -*-

;;; Commentary:

;; Runs the same checks as the Makefile without make or a Unix shell, so
;; they work under native Windows Emacs too:
;;
;;   emacs -Q --batch -l test/run-tests.el            ; ERT suite
;;   emacs -Q --batch -l test/run-tests.el compile    ; warnings are errors
;;   emacs -Q --batch -l test/run-tests.el checkdoc
;;   emacs -Q --batch -l test/run-tests.el all        ; all three
;;
;; Exits non-zero on any failure.

;;; Code:

(require 'cl-lib)

(defconst fastq-run--test-dir
  (file-name-directory (or load-file-name buffer-file-name))
  "The test directory.")

(defconst fastq-run--root
  (file-name-directory (directory-file-name fastq-run--test-dir))
  "The repository root.")

(defun fastq-run--sources ()
  "Return the package source files (absolute paths)."
  (directory-files fastq-run--root t "\\`fastq-.*\\.el\\'"))

(defun fastq-run--setup ()
  "Put the package and tests on `load-path'; never load stale .elc files."
  (setq load-prefer-newer t)
  (push fastq-run--root load-path)
  (push fastq-run--test-dir load-path))

(defun fastq-run-compile ()
  "Byte-compile every source with warnings as errors; return non-nil if ok."
  (require 'bytecomp)
  (let ((byte-compile-error-on-warn t)
        (ok t))
    (dolist (f (fastq-run--sources))
      (let ((elc (byte-compile-dest-file f)))
        (unwind-protect
            (unless (condition-case err (byte-compile-file f)
                      (error (message "%s: %s" f (error-message-string err)) nil))
              (message "byte-compile FAILED: %s" f)
              (setq ok nil))
          (when (file-exists-p elc) (delete-file elc)))))
    ok))

(defun fastq-run-checkdoc ()
  "Run checkdoc over every source; return non-nil if it reports nothing."
  (require 'checkdoc)
  (let ((checkdoc-diagnostic-buffer "*fastq-checkdoc*"))
    (dolist (f (fastq-run--sources))
      (with-current-buffer (find-file-noselect f)
        (checkdoc-current-buffer t)))
    (let ((b (get-buffer checkdoc-diagnostic-buffer)))
      (if (and b (with-current-buffer b
                   (goto-char (point-min))
                   (re-search-forward "^[^*\n\f]" nil t)))
          (progn (princ (with-current-buffer b (buffer-string))) nil)
        t))))

(defun fastq-run-tests ()
  "Load every test file and run the ERT suite; exit with its status."
  (dolist (f (directory-files fastq-run--test-dir t "-test\\.el\\'"))
    (load f nil t))
  (ert-run-tests-batch-and-exit))

(defun fastq-run-main ()
  "Dispatch on the remaining command-line arguments."
  (fastq-run--setup)
  (let ((what (or (car command-line-args-left) "test")))
    (setq command-line-args-left nil)
    (pcase what
      ("compile" (kill-emacs (if (fastq-run-compile) 0 1)))
      ("checkdoc" (kill-emacs (if (fastq-run-checkdoc) 0 1)))
      ("all" (unless (and (fastq-run-compile) (fastq-run-checkdoc))
               (kill-emacs 1))
       (fastq-run-tests))
      ("test" (fastq-run-tests))
      (_ (message "Usage: emacs -Q --batch -l test/run-tests.el [test|compile|checkdoc|all]")
         (kill-emacs 2)))))

(when noninteractive
  (fastq-run-main))

(provide 'run-tests)
;;; run-tests.el ends here
