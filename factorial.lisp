(defun fact (n)
  "simple implementation"
  (if (<= n 1) 1
      (* n (fact (1- n)))))

(declaim (ftype (function (fixnum) (values fixnum &optional)) fact-fn))
(defun fact-fn (n)
  "factorial calclation in type:fixnum"
  (declare (type fixnum n))
  (if (<= n 1) 1
      (* n (fact-fn (1- n)))))

(defun fact-tro (n)
  "tail recusion"
  (labels ((sub-fact (n acc)
	     (if (<= n 1) acc
		 (sub-fact (1- n) (* n acc)))))
    (sub-fact n 1)))

(defun fact-do (n)
  "not rucusive"
  (let ((result 1))
    (do ((i 1 (1+ i)))
	((> i n) result)
      (setq result (* result (1+ i))))))
