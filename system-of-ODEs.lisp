;;; system of first-order ordinary differential equations solvers

(setq *read-default-float-format* 'double-float)

;; make "graph-cl/" in current directory
(let ((dir (merge-pathnames #P"graph-cl/" (truename "./"))))
  (ensure-directories-exist dir))

(defun plot-first-order-odes (method list-funcs init-values interval-step num-of-steps &optional (marker-of-dat-file 0))
  (with-open-file (fp (format nil "~A~A~A-~A.dat"
			      (directory-namestring (truename "./"))
			      "graph-cl/"
			      "ODEs"
			      marker-of-dat-file)
		      :direction :output
		      :if-exists :supersede
		      :if-does-not-exist :create)
    (let ((curr-values init-values))
      (dotimes (i num-of-steps)
	(dolist (value curr-values) (format fp "~A~T" value))
	(format fp "~%")
	(setf curr-values
	      (funcall method list-funcs curr-values interval-step))))))

;; multi-euler method
(defun meuler-method (list-funcs curr-values interval-step)
  (let ((result-rev nil))
    (dolist (func list-funcs (mapcar #'+ curr-values (reverse result-rev)))
      (setf result-rev (cons (* interval-step (apply func curr-values)) result-rev)))))

;; multi-classical Runge--Kutta method
(defun mcrk-method (list-funcs curr-values interval-step)
  (let ((k1 nil)
	(k2 nil)
	(k3 nil)
	(k4 nil))
    (dolist (func list-funcs)		;k1
      (setf k1 (cons (* interval-step (apply func curr-values))	k1)))
    (setf k1 (reverse k1))
    (dolist (func list-funcs)		;k2
      (setf k2 (cons (* interval-step (apply func (mapcar (lambda (x y) (+ x (/ y 2.0))) curr-values k1))) k2)))
    (setf k2 (reverse k2))
    (dolist (func list-funcs)		;k3
      (setf k3 (cons (* interval-step (apply func (mapcar (lambda (x y) (+ x (/ y 2.0))) curr-values k2))) k3)))
    (setf k3 (reverse k3))
    (dolist (func list-funcs)		;k4
      (setf k4 (cons (* interval-step (apply func (mapcar #'+ curr-values k3))) k4)))
    (setf k4 (reverse k4))
    (mapcar (lambda (curr a b c d) (+ curr (/ (+ a (* 2.0 b) (* 2.0 c) d) 6.0)))
	    curr-values k1 k2 k3 k4)))

;; some functions for mrk-method
(defun multi-apply (funcs args)
  "(f1, ..., fn): (x1, ..., xn) |-> (f1(x1, ..., xn), ..., fn(x1, x2, ..., xn))"
  (let ((tmp '()))
    (dolist (func funcs (reverse tmp)) (push (apply func args) tmp))))

(defun mult-scalar (scalar args)
  "a: (x1, ..., xn) |-> (a*x1, ..., a*xn)"
  (mapcar (lambda (x) (* scalar x)) args))

(defun sum-mscalar (scalars lists)
  "(a1, ..., as): ((x11, ..., x1n), .., (xs1, ..., xsn)) |-> (a1*x11+...+as*xs1, ..., a1*x1n+...+as*xsn)"
  (reduce (lambda (x y) (mapcar #'+ x y))
	  (loop for scalar in scalars
	     for list in lists
	     collect (mult-scalar scalar list))
	  :initial-value (mult-scalar 0 (car lists))))

;; multi-Runge--Kutta method
(defun mrk-method (funcs curr-values interval-step &key (butcher-a '((1/2) (0 1/2) (0 0 1))) (butcher-b '(1/6 1/3 1/3 1/6)))
  (let ((k-list '()))
    (setf k-list (append k-list		;k1
			 (list (mult-scalar interval-step
					    (multi-apply funcs curr-values)))))
    (dolist (stage butcher-a)		;k2--kn
      (setf k-list (append k-list
			   (list (mult-scalar interval-step
					      (multi-apply funcs
							   (mapcar #'+ curr-values (sum-mscalar stage k-list))))))))
    (mapcar #'+ curr-values (sum-mscalar butcher-b k-list))))

;; Butcher's tableau
(defparameter *rk4-a-4s4o* '((1/2) (0 1/2) (0 0 1)))
(defparameter *rk4-b-4s4o* '(1/6 1/3 1/3 1/6))

(defparameter *shanks-a-9s7o* '((2/9) (1/12 1/4) (1/8 0 3/8) (23/216 0 7/72 -1/27) (-4136/729 0 -4528/243 5264/729 1456/81) (8087/11664 0 484/243 -518/729 -658/351 7/624) (-1217/2160 0 -145/72 8342/6615 361/195 3033/50960 117/490) (259/2768 0 -84/173 -14/173 6210/2249 -99873/251888 -29160/15743 2160/2249)))
(defparameter *shanks-b-9s7o* '(173/3360 0 0 1846/5145 27/91 -19683/713440 -19683/713440 27/91 173/3360))


;;; some examples
(defun example-crk-method-pendulum ()
  (labels ((pendulum-func (k)
	     (list (lambda (tau angle vel) (+ (* 0 tau angle vel) 1))
		   (lambda (tau angle vel) (+ (* 0 tau angle) vel))
		   (lambda (tau angle vel) (- (* 0 tau vel) (* k (sin angle)))))))
    (plot-first-order-odes #'mcrk-method (pendulum-func 1.0) '(0.0 -1.0 1.75) 5.0e-1 200 "cRK-pendulum")))

(defun example-rk-method-pendulum ()
  (labels ((pendulum-func (k)
	     (list (lambda (tau angle vel) (+ (* 0 tau angle vel) 1))
		   (lambda (tau angle vel) (+ (* 0 tau angle) vel))
		   (lambda (tau angle vel) (- (* 0 tau vel) (* k (sin angle))))))
	   (mrk (funcs curr-values interval-step)
	     (mrk-method funcs curr-values interval-step
			 :butcher-a *shanks-a-9s7d*
			 :butcher-b *shanks-b-9s7d*)))
    (plot-first-order-odes #'mrk (pendulum-func 1.0) '(0.0 -1.0 1.75) 5.0e-1 200 "RK-pendulum")))

(defun example-crk-method-mathieu ()
  "d^2x/dt^2 = b - alpha*x^n - (d + e*cos(t))*x"
  (labels ((mathieu-func (delta epsilon bias alpha n)
	     (list (lambda (tau pos vel) (+ (* 0 tau pos vel) 1))
		   (lambda (tau pos vel) (+ (* 0 tau pos) vel))
		   (lambda (tau pos vel) (- (+ (* 0 vel) bias) (* alpha (expt pos n))
					    (* pos (+ delta (* epsilon (cos tau)))))))))
    (plot-first-order-odes #'mcrk-method (mathieu-func 1.5 10.0 0.5 1.0 3) '(0.0 1.0 0.0) 5.0e-2 4000 "cRK-mathieu")))

(defun example-rk-method-mathieu ()
  "d^2x/dt^2 = b - alpha*x^n - (d + e*cos(t))*x"
  (labels ((mathieu-func (delta epsilon bias alpha n)
	     (list (lambda (tau pos vel) (+ (* 0 tau pos vel) 1))
		   (lambda (tau pos vel) (+ (* 0 tau pos) vel))
		   (lambda (tau pos vel) (- (+ (* 0 vel) bias) (* alpha (expt pos n))
					    (* pos (+ delta (* epsilon (cos tau))))))))
	   (mrk (funcs curr-values interval-step)
	     (mrk-method funcs curr-values interval-step
			 :butcher-a *shanks-a-9s7d*
			 :butcher-b *shanks-b-9s7d*)))
    (plot-first-order-odes #'mrk (mathieu-func 1.5 10.0 0.5 1.0 3) '(0.0 1.0 0.0) 5.0e-2 4000 "RK-mathieu")))
