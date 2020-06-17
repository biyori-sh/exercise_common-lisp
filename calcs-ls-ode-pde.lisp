(setq *read-default-float-format* 'double-float)
;; (ql:quickload :lparallel)
;; (setf lparallel:*kernel* (lparallel:make-kernel 2))

;;; line searches
(defun golden-divide-method (init end func &optional (eps 1.0e-6))
  "search for a minimum of downward-convex function in a range from init to end"
  (let* ((left-edge init)
	 (right-edge end)
	 (ratio-golden (/ (- (sqrt 5.0e0) 1.0e0) 2.0e0))
	 (left-inner (+ (* ratio-golden left-edge)
			(* (- 1 ratio-golden) right-edge)))
	 (right-inner (+ (* (- 1 ratio-golden) left-edge)
			 (* ratio-golden right-edge)))
	 (value-left (funcall func left-inner))
	 (value-right (funcall func right-inner))
	 (count-loop 0)
	 (flag nil)
	 (result 0.0))
    (loop (when (< (- right-edge left-edge) eps) (return nil))
       (setq value-left (funcall func left-inner)
	     value-right (funcall func right-inner))
       (incf count-loop)
       (cond ((< value-left value-right)
	      (setq right-edge right-inner
		    right-inner left-inner
		    left-inner (+ (* ratio-golden left-edge)
				  (* (- 1 ratio-golden) right-edge))))
	     ((>= value-left value-right)
	      (setq left-edge left-inner
		    left-inner right-inner
		    right-inner (+ (* (- 1 ratio-golden) left-edge)
				   (* ratio-golden right-edge))))))
    (cond ((< (abs (- left-edge init)) (* 10 eps)) (setq result init
						  flag -1))
	  ((< (abs (- end right-edge)) (* 10 eps)) (setq result end
						  flag 1))
	  (t (setq result (/ (+ left-edge right-edge) 2.0)
		   flag nil)))
    (values result count-loop flag)))


(defun bisection-method (init end func &optional (eps 1.0e-6))
  "find root by using bisection method in a range from init to end"
  (let* ((edge-l init)
	 (edge-r end)
	 (middle (/ (+ init end) 2.0))
	 (value-l (funcall func init))
	 (value-r (funcall func end))
	 (value-m (funcall func middle)))
    (when (minusp (* value-l value-r))
      (loop (when (< (- edge-r edge-l) eps)
	      (return (/ (+ edge-l edge-r) 2.0)))
	 (if (minusp (* value-m value-r))
	     (progn (setf edge-l middle
			  middle (/ (+ edge-l edge-r) 2.0)
			  value-l (funcall func edge-l)
			  value-m (funcall func middle)))
	     (progn (setf edge-r middle
			  middle (/ (+ edge-l edge-r) 2.0)
			  value-m (funcall func middle)
			  value-r (funcall func edge-r))))))))

(defun newtons-method-primitive (init-appro func &optional (eps 1.0e-6))
  "primitive implementation of Newton's method starting from init-appro"
  (let ((curr init-appro)
	(next 0.0))
    (labels ((deriv-func (point func)
	       (/ (- (funcall func (+ point eps))
		     (funcall func (- point eps)))
		  (* 2 eps)))
	     (next-appro (curr-appro func)
	       (- curr-appro (/ (funcall func curr-appro)
				(deriv-func curr-appro func)))))
      (setf next (next-appro curr func))
      (loop (when (< (abs (- curr next)) eps) (return next))
	 (setf curr next
	       next (next-appro curr func))))))


;;; first-order ordinary differential equation solvers
;; make "graph-cl/" in current directory
(let ((dir (merge-pathnames #P"graph-cl/" (truename "./"))))
  (ensure-directories-exist dir))

(defun plot-first-order-ode (method func-of-time-value init-value interval-time num-of-steps &optional (marker-of-dat-file 0))
  "solve first-order ordinary differential equation, dx/dt = f(t, x), using designated method and plot data to file"
  (with-open-file (fp (format nil "~A~A~A-~A-~A.dat"
			      (directory-namestring (truename "./"))
			      "graph-cl/"
			      "ODE-FO"
			      method
			      marker-of-dat-file)
		      :direction :output
		      :if-exists :supersede
		      :if-does-not-exist :create)
    (let ((curr-time 0.0)
	  (curr-value init-value))
      (dotimes (i num-of-steps)
	(format fp "~A~T~A~%"  curr-time curr-value)
	(setf (values curr-time curr-value)
	      (funcall method func-of-time-value curr-time curr-value interval-time))))))

(defun euler-method (func-of-time-value curr-time curr-value interval-time)
  "calculate next (time, value) using Euler method for first-order ordinary differential equation, dx/dt = f(t, x)"
  (values (+ curr-time interval-time)
	  (+ curr-value (* interval-time (funcall func-of-time-value
						  curr-time
						  curr-value)))))

(defun implicit-euler-method (func-of-time-value curr-time curr-value interval-time)
  "calculate next (time, value) using implicit Euler method for first-order ordinary differential equation; be implemented using newtons-method-primitive (central difference) and can occur zero-division in some situations (insufficient accuracy,  too large or small derivative, and ...)"
  (let ((next-time (+ interval-time curr-time)))
    (values next-time
	    (newtons-method-primitive
	     curr-value
	     (lambda (next-value)
	       (- next-value curr-value (* interval-time
					   (funcall func-of-time-value
						    next-time
						    next-value))))))))

(defun runge-kutta-method (func-of-time-value curr-time curr-value interval-time)
  "calculate next (time, value) using fourth-order Runge--Kutta method with four-stages for first-order ordinary differential equation, dx/dt = f(t, x)"
  (let* ((k_1 (* interval-time (funcall func-of-time-value
					curr-time
					curr-value)))
	 (k_2 (* interval-time (funcall func-of-time-value
					(+ curr-time (/ interval-time 2.0))
					(+ curr-value (/ k_1 2.0)))))
	 (k_3 (* interval-time (funcall func-of-time-value
					(+ curr-time (/ interval-time 2.0))
					(+ curr-value (/ k_2 2.0)))))
	 (k_4 (* interval-time (funcall func-of-time-value
					(+ curr-time interval-time)
					(+ curr-value k_3)))))
    (values (+ curr-time interval-time)
	    (+ curr-value (/ (+ k_1 (* 2 k_2) (* 2 k_3) k_4) 6.0)))))

;; some examples
(defun example-euler-method ()
  (labels ((test-func ()
	     (lambda (val-t val-x) (* 500 val-x val-x (- 1 val-x)))))
    (plot-first-order-ode 'euler-method (test-func) (float 1/100) (float 1/200) 200 200)
    (plot-first-order-ode 'euler-method (test-func) (float 1/100) (float 1/400) 400 400)))

(defun example-implicit-euler-method ()
  (labels ((test-func ()
	     (lambda (val-t val-x) (* 500 val-x val-x (- 1 val-x)))))
    (plot-first-order-ode 'implicit-euler-method (test-func) (float 1/100) (float 1/200) 200 200)
    (plot-first-order-ode 'implicit-euler-method (test-func) (float 1/100) (float 1/400) 400 400)))

(defun example-runge-kutta-method ()
  (labels ((test-func ()
	     (lambda (val-t val-x) (* 500 val-x val-x (- 1 val-x)))))
    (plot-first-order-ode 'runge-kutta-method (test-func) (float 1/100) (float 1/200) 200 200)
    (plot-first-order-ode 'runge-kutta-method (test-func) (float 1/100) (float 1/400) 400 400)))
