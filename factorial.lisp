;; (declaim (ftype (function (fixnum) (values fixnum &optional)) fact))
(defun fact-r (n)
  (declare (type fixnum n))
  (if (<= n 1) 1
      (* n (fact-r (1- n)))))

(defun fact-tro (n)
  (declare (type fixnum n))
  (labels ((sub-fact (n acc)
	     (if (<= n 1) acc
		 (sub-fact (1- n) (* n acc)))))
    (sub-fact n 1)))
