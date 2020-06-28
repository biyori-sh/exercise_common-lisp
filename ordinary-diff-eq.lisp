;;; first-order ordinary differential equation solvers for one dimension

;; for implicit euler method (newtons-method-primitive)
(load "optimization.lisp")

;; make "graph-cl/" in current directory
(let ((dir (merge-pathnames #P"graph-cl/" (truename "./"))))
  (ensure-directories-exist dir))

(defun plot-first-order-ode (method func-of-time-position init-position interval-time num-of-steps &optional (marker-of-dat-file 0))
  "solve first-order ordinary differential equation, dx/dt = f(t, x), using designated method and plot data to file"
  (with-open-file (fp (format nil "~A~A~A-~A-~A.dat"
			      (directory-namestring (truename "./"))
			      "graph-cl/"
			      "first-order-ODE"
			      method
			      marker-of-dat-file)
		      :direction :output
		      :if-exists :supersede
		      :if-does-not-exist :create)
    (let ((curr-time 0.0)
	  (curr-position init-position))
      (dotimes (i num-of-steps)
	(format fp "~A~T~A~%"  curr-time curr-position)
	(setf (values curr-time curr-position)
	      (funcall method func-of-time-position curr-time curr-position interval-time))))))

;; explicit euler method
(defun euler-method (func-of-time-position curr-time curr-position interval-time)
  "calculate next (time, value) using Euler method for first-order ordinary differential equation, dx/dt = f(t, x)"
  (values (+ curr-time interval-time)
	  (+ curr-position (* interval-time (funcall func-of-time-position
						     curr-time
						     curr-position)))))
;; implicit euler method
(defun implicit-euler-method (func-of-time-position curr-time curr-position interval-time)
  "calculate next (time, value) using implicit Euler method for first-order ordinary differential equation; be implemented using newtons-method-primitive (central difference) and can occur zero-division in some situations (insufficient accuracy,  too large or small derivative, and ...)"
  (let ((next-time (+ interval-time curr-time)))
    (values next-time
	    (newtons-method-primitive
	     curr-position
	     (lambda (next-position)
	       (- next-position curr-position (* interval-time
						 (funcall func-of-time-position
							  next-time
							  next-position))))))))

;; classical Runge--Kutta method
(defun runge-kutta-method (func-of-time-position curr-time curr-position interval-time)
  "calculate next (time, value) using fourth-order Runge--Kutta method with four-stages for first-order ordinary differential equation, dx/dt = f(t, x)"
  (let* ((k_1 (* interval-time (funcall func-of-time-position
					curr-time
					curr-position)))
	 (k_2 (* interval-time (funcall func-of-time-position
					(+ curr-time (/ interval-time 2.0))
					(+ curr-position (/ k_1 2.0)))))
	 (k_3 (* interval-time (funcall func-of-time-position
					(+ curr-time (/ interval-time 2.0))
					(+ curr-position (/ k_2 2.0)))))
	 (k_4 (* interval-time (funcall func-of-time-position
					(+ curr-time interval-time)
					(+ curr-position k_3)))))
    (values (+ curr-time interval-time)
	    (+ curr-position (/ (+ k_1 (* 2 k_2) (* 2 k_3) k_4) 6.0)))))


;; some examples
(defun example-euler-method ()
  (labels ((test-func ()
	     (lambda (val-t val-x) (* 500 val-x val-x (- 1 val-x (* 0 val-t))))))
    (plot-first-order-ode 'euler-method (test-func) (float 1/100) (float 1/200) 200 200)
    (plot-first-order-ode 'euler-method (test-func) (float 1/100) (float 1/400) 400 400)))

(defun example-implicit-euler-method ()
  (labels ((test-func ()
	     (lambda (val-t val-x) (* 500 val-x val-x (- 1 val-x (* 0 val-t))))))
    (plot-first-order-ode 'implicit-euler-method (test-func) (float 1/100) (float 1/200) 200 200)
    (plot-first-order-ode 'implicit-euler-method (test-func) (float 1/100) (float 1/400) 400 400)))

(defun example-runge-kutta-method ()
  (labels ((test-func ()
	     (lambda (val-t val-x) (* 500 val-x val-x (- 1 val-x (* 0 val-t))))))
    (plot-first-order-ode 'runge-kutta-method (test-func) (float 1/100) (float 1/200) 200 200)
    (plot-first-order-ode 'runge-kutta-method (test-func) (float 1/100) (float 1/400) 400 400)))
