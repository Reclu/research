;; (require 'tex) ; pour être sur d'avoir le mode TeX
;; (require 'latex) ; pour être sur d'avoir le mode LaTeX
(require 'reftex)

;; no startup msg  
(setq inhibit-startup-message t)        ; Disable startup message

(setq package-archives '(("gnu" . "http://elpa.gnu.org/packages/")  ("melpa" . "http://melpa.milkbox.net/packages/")))

(setq package-enable-at-startup nil)
  (package-initialize)

; Turn on RefTeX for AUCTeX, http://www.gnu.org/s/auctex/manual/reftex/reftex_5.html
(add-hook 'LaTeX-mode-hook 'turn-on-reftex)
; Make RefTeX interact with AUCTeX, http://www.gnu.org/s/auctex/manual/reftex/AUCTeX_002dRefTeX-Interface.html
(setq reftex-plug-into-AUCTeX t)

;; LaTeX with tikz external
(eval-after-load "tex"
  '(add-to-list 'TeX-command-list 
		'("LaTeX-Ext" "pdflatex -shell-escape %s" TeX-run-TeX nil t :help "Run LaTeX with tikz external")))


;; (add-hook 'LaTeX-mode-hook 'TeX-PDF-mode) ;turn on pdf-mode.  AUCTeX
(add-hook 'LaTeX-mode-hook 'LaTeX-PDF-mode) ;turn on pdf-mode.  AUCTeX
(add-hook 'LaTeX-mode-hook 'LaTeX-math-mode) ;turn on math-mode by default
(add-hook 'LaTeX-mode-hook 'reftex-mode) ;turn on REFTeX mode by default
(add-hook 'LaTeX-mode-hook 'flyspell-mode) ;turn on flyspell mode by default

(add-hook 'LaTeX-mode-hook
      '(lambda ()
         (TeX-add-symbols '("eqref" TeX-arg-ref (ignore)))))
(add-hook 'LaTeX-mode-hook
      '(lambda ()
         (TeX-add-symbols '("subref" TeX-arg-ref (ignore)))))

(setq reftex-label-alist
  '((nil ?s nil nil nil ("Kapitel" "Kap." "Abschnitt" "Teil"))
    (nil ?e nil nil nil ("Gleichung" "Gl."))
    (nil ?t nil nil nil ("Tabelle"))
    (nil ?f nil nil nil ("Figur" "Abbildung" "Abb."))
    (nil ?n nil nil nil ("Anmerkung" "Anm."))
    (nil ?i nil nil nil ("Punkt"))))

(setq reftex-label-alist
      '(("example"   ?a "ex:"  "~\\ref{%s}" nil ("example"   "ex.") -2)
	("remark"   ?r "rq:"  "~\\ref{%s}" nil ("remark"   "rk.") -3)
          ("theorem" ?h "th:" "~\\ref{%s}" t   ("theorem" "th.") -4)))

(add-hook 'LaTeX-mode-hook
        (lambda ()
          (LaTeX-add-environments
	   '("example" LaTeX-env-label)
	   '("remark" LaTeX-env-label)
            '("theorem" LaTeX-env-label))))


;; (setq reftex-label-alist
;;   '(("equation"   ?e "eq:"  "(~\\ref{%s})" nil ("Gleichung"   "Gl.")))


;(setq reftex-label-alist '((nil ?e nil "~\\eqref{%s}" nil nil)))
;(setq reftex-label-alist '(AMSTeX))
     
(require 'auto-complete-config)
(ac-config-default)
(setq ac-auto-show-menu nil)
(ac-set-trigger-key "TAB")
(add-to-list 'ac-modes 'latex-mode)
 
(require 'auto-complete-auctex) ;
