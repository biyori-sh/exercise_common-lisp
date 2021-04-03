;; Why You Cannot (Yet) Write an “Interval Arithmetic” Library in Common Lisp
;; . . . or: Hammering Some Sense into :ieee-floating-point
;; Marco Antoniotti
;; https://arxiv.org/pdf/2003.03831.pdf

(defparameter *default-rounding-mode*
  (second (member :rounding-mode (sb-int:get-floating-point-modes))))


(defmacro round-to (mode-keyword &body body)
  "Evaluate BODY in the rounding mode, :POSITIVE-INFINITY or :NEGATIVE-INFINITY or :NEAREST. After that, the rounding-mode returns to *default-rounding-mode*."
  (labels ((rounding-modep (mode)
             (declare (inline rounding-modep)
                      (keyword mode))
             (find mode (list :NEAREST :POSITIVE-INFINITY :NEGATIVE-INFINITY))))
    (when (rounding-modep mode-keyword)
      `(prog2 (sb-int:set-floating-point-modes :ROUNDING-MODE ,mode-keyword)
           (progn ,@body)
         (sb-int:set-floating-point-modes :ROUNDING-MODE *default-rounding-mode*)))))

(defmacro round-positive (&body body)
  "Evaluate BODY in the rounding mode: POSITIVE-INFINITY. After that, the rounding-mode returns to *default-rounding-mode*."
  `(prog2 (sb-int:set-floating-point-modes :ROUNDING-MODE :POSITIVE-INFINITY)
       (progn ,@body)
     (sb-int:set-floating-point-modes :ROUNDING-MODE *default-rounding-mode*)))

(defmacro round-negative (&body body)
  "Evaluate BODY in the rounding mode: NEGATIVE-INFINITY. After that, the rounding-mode returns to *default-rounding-mode*."
  `(prog2 (sb-int:set-floating-point-modes :ROUNDING-MODE :NEGATIVE-INFINITY)
       (progn ,@body)
     (sb-int:set-floating-point-modes :ROUNDING-MODE *default-rounding-mode*)))

(defmacro round-nearest (&body body)
  "Evaluate BODY in the rounding mode: NEAREST. After that, the rounding-mode returns to *default-rounding-mode*."
  `(prog2 (sb-int:set-floating-point-modes :ROUNDING-MODE :NEAREST)
       (progn ,@body)
     (sb-int:set-floating-point-modes :ROUNDING-MODE *default-rounding-mode*)))


(defstruct ([] (:constructor [] (low high)))
  "Structure type for the interval arithmetic."
  (low 0.0 :type real)
  (high 0.0 :type real))

(defmethod diameteri ((i []))
  (- ([]-high i) ([]-low i)))

(defmethod radiusi ((i []))
  (/ (diameteri i) 2.0))

(defmethod centeri ((i []))
  (/ (- ([]-high i) ([]-low i)) 2.0))

(defmethod pointip ((i []))
  (= ([]-high i) ([]-low i)))

(defmethod point->interval ((point real))
  ([] point point))

(defmethod elementip ((point number) (i []))
  (and (<= ([]-low i) point) (<= point ([]-high i))))

(defmethod subintervalip ((subi []) (superi []))
  (and (<= ([]-low superi) ([]-low subi)) (<= ([]-high subi) ([]-high superi))))

(defmethod zeroip ((i []))
  (elementip 0 i))

(defmethod addi ((i1 []) (i2 []))
  ([] (round-negative (+ ([]-low i1) ([]-low i2)))
      (round-positive (+ ([]-high i1) ([]-high i2)))))

(defmethod subi ((i1 []) (i2 []))
  ([] (round-negative (- ([]-low i1) ([]-high i2)))
      (round-positive (- ([]-high i1) ([]-low i2)))))

(defmethod muli ((i1 []) (i2 []))
  ([] (round-negative
        (min (* ([]-high i1) ([]-high i2))
             (* ([]-low i1) ([]-high i2))
             (* ([]-high i1) ([]-low i2))
             (* ([]-low i1) ([]-low i2))))
      (round-positive
        (max (* ([]-high i1) ([]-high i2))
             (* ([]-low i1) ([]-high i2))
             (* ([]-high i1) ([]-low i2))
             (* ([]-low i1) ([]-low i2))))))

(defmethod divi ((i1 []) (i2 []))
  (if (zeroip i2)
      (error "Not defined for the divisor including zero: ~a~%" i2)
      (muli i1 ([] (/ 1.0 ([]-high i2)) (/ 1.0 ([]-low i2))))))

;; (sb-int:get-floating-point-modes)
;; (sb-int:set-floating-point-modes :ROUNDING-MODE :POSITIVE-INFINITY)
;; (sb-int:set-floating-point-modes :ROUNDING-MODE :NEGATIVE-INFINITY)
;; (sb-int:set-floating-point-modes :ROUNDING-MODE :NEAREST)
