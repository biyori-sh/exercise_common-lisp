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
