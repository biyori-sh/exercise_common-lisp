;; Why You Cannot (Yet) Write an “Interval Arithmetic” Library in Common Lisp
;; . . . or: Hammering Some Sense into :ieee-floating-point
;; Marco Antoniotti
;; https://arxiv.org/pdf/2003.03831.pdf

;; only fo SBCL
(defparameter *default-rounding-mode*
  (getf (sb-int:get-floating-point-modes) :ROUNDING-MODE)
  "a keyword of the default rounding mode")

;; Macro for rounding mode
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


;; type
(deftype positive-real ()
  '(real (0) *))

(deftype negative-real ()
  '(real * (0)))


;; Structure type of Interval
;; check large/small relation between low and high, required functions?
;; more detailed? positive interval, interval including zero and negative interval
(defstruct ([] (:constructor [] (low high)))
  "Structure type for the interval arithmetic."
  (low 0.0 :type real)
  (high 0.0 :type real))


;; predicates
(defmethod pointp ((i []))
  (= ([]-high i) ([]-low i)))

(defmethod wd[]p ((i []))
  "Is this a well-defined interval?"
  (<= ([]-low i) ([]-high i)))

(defmethod elementp ((point number) (i []))
  (and (wd[]p i)
       (<= ([]-low i) point) (<= point ([]-high i))))

(defmethod sub[]p ((subi []) (superi []))
  (and (wd[]p subi) (wd[]p superi)
       (<= ([]-low superi) ([]-low subi)) (<= ([]-high subi) ([]-high superi))))

(defmethod zero[]p ((i []))
  (and (wd[]p i) (elementp 0 i)))

(defmethod positive[]p ((i []))
  (and (wd[]p i) (< 0 ([]-low i))))

(defmethod negative[]p ((i []))
  (and (wd[]p i) (< ([]-high i) 0)))

;; Unary operators
(defmethod diameter ((i []))
  (- ([]-high i) ([]-low i)))

(defmethod radius ((i []))
  (/ (diameter i) 2.0))

(defmethod center ((i []))
  (/ (- ([]-high i) ([]-low i)) 2.0))

(defmethod point->[] ((point real))
  ([] point point))

(defmethod abs[]-min ((i []))
  (if (zero[]p i) 0
      (min (abs ([]-low i)) (abs ([]-high i)))))

(defmethod abs[]-max ((i []))
  (max (abs ([]-low i)) (abs ([]-high i))))

(defmethod abs[] ((i []))
  ([] (round-negative (abs[]-min i))
      (round-positive (abs[]-max i))))

(defmethod abs[] ((i real))
  (point->[] (abs i)))

;; Binary operators
;; Addittion
(defmethod +[] ((i1 []) (i2 []))
  ([] (round-negative (+ ([]-low i1) ([]-low i2)))
      (round-positive (+ ([]-high i1) ([]-high i2)))))

(defmethod +[] ((i1 real) (i2 []))
  ([] (round-negative (+ i1 ([]-low i2)))
      (round-positive (+ i1 ([]-high i2)))))

(defmethod +[] ((i1 []) (i2 real))
  (add[] i2 i1))

(defmethod +[] ((i1 real) (i2 real))
  ([] (round-negative (+ i1 i2))
      (round-positive (+ i1 i2))))

;; Subtraction
(defmethod -[] ((i1 []) (i2 []))
  (+[] i1 ([] (- ([]-low i2)) (- ([]-high i2)))))

(defmethod -[] ((i1 real) (i2 []))
  (+[] (point->[] i1) ([] (- ([]-low i2)) (- ([]-high i2)))))

(defmethod -[] ((i1 []) (i2 real))
  (+[] i1 (point->[] (- i2))))

(defmethod -[] ((i1 real) (i2 real))
  (+[] (point->[] i1) (point->[] (- i2))))

;; Multiplication
;; (defmethod *[] ((i1 []) (i2 []))
;;   ([] (round-negative
;;         (cond ((and (positive[]p i1) (positive[]p i2))
;;                (* ([]-low i1) ([]-low i2)))
;;               ((and (negative[]p i1) (negative[]p i2))
;;                (* ([]-high i1) ([]-high i2)))
;;               (t (min (* ([]-low i1) ([]-high i2))
;;                       (* ([]-high i1) ([]-low i2))))))
;;       (round-positive
;;         (cond ((and (positive[]p i1) (negative[]p i2))
;;                (* ([]-low i1) ([]-high i2)))
;;               ((and (negative[]p i1) (positive[]p i2))
;;                (* ([]-high i1) ([]-low i2)))
;;               (t (max (* ([]-high i1) ([]-high i2))
;;                       (* ([]-low i1) ([]-low i2))))))))

(defmethod *[] ((i1 []) (i2 []))
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

(defmethod *[] ((i1 real) (i2 []))
  ([] (round-negative
        (min (* i1 ([]-high i2))
             (* i1 ([]-low i2))))
      (round-positive
        (max (* i1 ([]-high i2))
             (* i1 ([]-low i2))))))

(defmethod *[] ((i1 []) (i2 real))
  (*[] i2 i1))

(defmethod *[] ((i1 real) (i2 real))
  ([] (round-negative (* i1 i2))
      (round-positive (* i1 i2))))

;; Division
(defmethod /[] ((i1 []) (i2 []))
  (if (zero[]p i2)
      (error "Not defined for the divisor including zero: ~a~%" i2)
      (*[] i1 ([] (/ 1.0 ([]-high i2)) (/ 1.0 ([]-low i2))))))

(defmethod /[] ((i1 real) (i2 []))
  (if (zero[]p i2)
      (error "Not defined for the divisor including zero: ~a~%" i2)
      (*[] (point->[] i1) ([] (/ 1.0 ([]-high i2)) (/ 1.0 ([]-low i2))))))

(defmethod /[] ((i1 []) (i2 real))
  (*[] i1 (point->[] (/ 1.0 i2))))

(defmethod /[] ((i1 real) (i2 real))
  (*[] (point->[] i1) (point->[] (/ 1.0 i2))))

;; (let ((result ([] 0.0 0.0)))
;;   (dotimes (i 100 result)
;;     (setf result (+[] 0.1 result))))
;; => #S([] :LOW 9.999982 :HIGH 10.000027)

;; (let ((result ([] 1.0 1.0)))
;;         (dotimes (i 100 result)
;;           (setf result (*[] ([] 0.95 0.95)  result))))
;; => #S([] :LOW 0.0059204977 :HIGH 0.00592055)

;; (sb-int:get-floating-point-modes)
;; (sb-int:set-floating-point-modes :ROUNDING-MODE :POSITIVE-INFINITY)
;; (sb-int:set-floating-point-modes :ROUNDING-MODE :NEGATIVE-INFINITY)
;; (sb-int:set-floating-point-modes :ROUNDING-MODE :NEAREST)
