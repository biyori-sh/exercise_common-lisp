(defun fib (n)
  "simple implementation"
  (cond ((<= n 0) 0)
	((= n 1) 1)
	(t (+ (fib (- n 1))
	      (fib (- n 2))))))

(declaim (inline fib-opt) (ftype (function (integer) integer) fib-opt))
(defun fib-opt (n)
  "practice to optimize"
  (if (<= n 1) n
      (+ (fib (- n 1)) (fib (- n 2)))))

(defun fib-tr (n)
  "tail recusion"
  (labels ((ftr-sub (cnter prev curr)
	     (cond ((= cnter 0) prev)
		   ((= cnter 1) curr)
		   (t (ftr-sub (1- cnter) curr (+ prev curr))))))
    (ftr-sub n 0 1)))

(defun fib-matrix (n)
  "fibonacci matrix"
  (labels ((product-matrix (ma mb)
	     (let ((result (make-array '(2 2) :initial-element 0)))
	       (dotimes (i 2)
		 (dotimes (j 2)
		   (let ((res (aref result i j)))
		     (dotimes (k 2)
		       (setf res (+ res (* (aref ma i k) (aref mb k j)))))
		     (setf (aref result i j) res))))
	       result)))
    (let ((fmatrix (make-array '(2 2) :initial-contents '((1 1) (1 0))))
	  (identity (make-array '(2 2) :initial-contents '((1 0) (0 1)))))
      (dotimes (i n)
	(setf identity (product-matrix fmatrix identity)))
      (values (aref identity 0 1) identity))))

(defun binary-lst (n)
  (mapcar #'digit-char-p (concatenate 'list (format nil "~b" n))))
