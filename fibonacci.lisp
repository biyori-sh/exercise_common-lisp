(defun fib (n)
  "simple implementation using recusion"
  (cond ((<= n 0) 0)
	((= n 1) 1)
	(t (+ (fib (- n 1))
	      (fib (- n 2))))))

(declaim (inline fib-opt) (ftype (function (integer) integer) fib-opt))
(defun fib-opt (n)
  "practice to optimize"
  (declare (type integer n)
	   (optimize (speed 3) (safety 0) (debug 0)))
  (the integer (if (<= n 1) n
                   (+ (fib-opt (- n 1)) (fib-opt (- n 2))))))

(defun fib-tr (n)
  "tail recusion"
  (labels ((ftr-sub (cnter prev curr)
	     (cond ((= cnter 0) prev)
		   ((= cnter 1) curr)
		   (t (ftr-sub (1- cnter) curr (+ prev curr))))))
    (ftr-sub n 0 1)))

(defun fib-dotimes (n)
  "not recusive"
  (let ((fib-prev 0)
	(fib-curr 1))
    (cond ((= n 0) fib-prev)
	  ((= n 1) fib-curr)
	  (t (dotimes (index (1- n) fib-curr)
	       (setf (values fib-prev fib-curr)
		     (values fib-curr (+ fib-prev fib-curr))))))))

(labels ((prod-mat (ma mb)
	   (let ((result (make-array '(2 2) :initial-element 0)))
	     (dotimes (i 2)
	       (dotimes (j 2)
		 (let ((res (aref result i j)))
		   (dotimes (k 2)
		     (setf res (+ res (* (aref ma i k) (aref mb k j)))))
		   (setf (aref result i j) res))))
	     result))
	 (make-binary-list (n)
	   (mapcar #'digit-char-p (concatenate 'list (format nil "~b" n)))))

  (defun fib-matrix (n)
    "calculate by using fibonacci matrix"
    (let ((mat-f (make-array '(2 2) :initial-contents '((1 1) (1 0))))
	  (mat-result (make-array '(2 2) :initial-contents '((1 0) (0 1)))))
      (dotimes (i n)
	(setf mat-result (prod-mat mat-f mat-result)))
      (values (aref mat-result 0 1) mat-result)))

  (defun fib-matrix-binary (n)
    "fib-matrix optimized by using binary list of power"
    (let ((lst-bin-rev (reverse (make-binary-list n)))
	  (mat-fib (make-array '(2 2) :initial-contents '((1 1) (1 0))))
	  (mat-result (make-array '(2 2) :initial-contents '((1 0) (0 1)))))
      (when (= 1 (car lst-bin-rev)) (prod-mat mat-fib mat-result))
      (dolist (index (cdr lst-bin-rev))
	(setf mat-fib (prod-mat mat-fib mat-fib))
	(when (= index 1)
	  (setf mat-result (prod-mat mat-fib mat-result))))
      (values (aref mat-result 0 1) mat-result))))


;;; generalized fibonacci
(defun fib-general (n list-inits)
  (let* ((num-preceding (length list-inits))
	 (list-fib (copy-list list-inits))
	 (len-list (length list-fib)))
    (loop (when (< n len-list)
	    (return (nth n list-fib)))
       (setf list-fib
	     (append list-fib
		     (list (reduce #'+ (last list-fib num-preceding)))))
       (incf len-list))))

(defun fib-general-rev (n list-inits)
  (labels ((sum-of-heads (k lst)
	     (let ((acc 0))
	       (dotimes (i k acc) (incf acc (nth i lst))))))
    (let* ((num-preceding (length list-inits))
	   (list-fib (reverse list-inits))
	   (len-list (length list-fib)))
      (loop (when (< n len-list)
	      (return (nth (- len-list n 1) list-fib)))
	 (push (sum-of-heads num-preceding list-fib) list-fib)
	 (incf len-list)))))

(declaim (inline fib-gen-trec)
	 (ftype (function (integer list) (integer)) fib-gen-trec))
(defun fib-gen-trec (n list-inits)
  (declare (type integer n)
	   (type list list-inits))
  (if (< n (length list-inits))
      (nth n list-inits)
      (fib-gen-trec (1- n) (append (cdr list-inits)
				  (list (reduce #'+ list-inits))))))
