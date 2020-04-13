;; (setq *read-default-float-format* 'double-float)

;;; line searches
(defun golden-divide-method (init end func &key (eps 1.0d-10))
  "searching for a minimum of downward-convex function in a range from init to end"
  (let* ((left-edge init)
	 (right-edge end)
	 (ratio-golden (/ (- (sqrt 5.0d0) 1.0d0) 2.0d0))
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


(defun bisection-method (init end func &key (eps 1.0d-10))
  "finding root by using bisection method in a range from init to end"
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


(defun newtons-method-primitive (init-appro func &key (eps 1.0d-10))
  "primitive implementation of newton's method starting from init-appro"
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


;;; ordinary differential equation solvers
(defun ode-euler-method (func-t-x init-value interval num-of-steps &optional (marker-of-dat-file 0))
  "solving ordinary differential equation, dx/dt = f(t, x), by using euler method"
  (with-open-file (fp (format nil "~A~A~A~A" (directory-namestring (truename ".")) "ode-euler-method-" marker-of-dat-file ".dat")
		      :direction :output
		      :if-exists :supersede
		      :if-does-not-exist :create)
    (let ((curr-time 0.0)
	  (curr-value init-value))
      (dotimes (i num-of-steps)
	(format fp "~a~t~a~%"  curr-time curr-value)
	(incf curr-time interval)
	(setf curr-value (+ curr-value (* interval (funcall func-t-x curr-time curr-value))))))))

(defun test-euler-method ()
  (labels ((test-func ()
	     (lambda (val-t val-x) (* 500 val-x val-x (- 1 val-x)))))
    (ode-euler-method (test-func) (float 1/100) (float 1/250) 250 250)
    (ode-euler-method (test-func) (float 1/100) (float 1/500) 500 500)))
