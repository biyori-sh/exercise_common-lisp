(defun fact (n)
  (declare (type fixnum n))
  (if (<= n 1) 1
       (* n (fact (1- n)))))
