(defun gen-point (&rest ranges-from-zero)
  "generate a random point of an n-dimensional space. each value is not negative."
  (labels ((gen-point-sub (coordinate ranges-from-zero)
             (if (null ranges-from-zero)
                 (reverse coordinate)
                 (gen-point-sub (cons (random (car ranges-from-zero)) coordinate)
                                (cdr ranges-from-zero)))))
    (gen-point-sub nil ranges-from-zero)))

(defun length-from-origin (point)
  "return the length between the point and the origin."
  (sqrt (reduce #'+ (mapcar (lambda (x) (* x x)) point))))

(defun monte-carlo-pi (num-of-trials)
  "calculate the pi by the monte carlo simulation; generate random points in a unit square and count the points in an unit quarter circle (return the ratio multiplied by four)."
  (let ((cnter 0))
    (dotimes (i num-of-trials)
      (when (<= (length-from-origin (gen-point 1.0 1.0)) 1.0)
        (incf cnter)))
    (* 4 (/ (* 1.0 cnter) num-of-trials))))

(defun parameters-pi (num-of-trials num-of-samplings)
  "run (monte-carlo-pi num-of-trials) by num-of-samplings times and return the average and variance."
  (let ((lst nil)
        (average 0.0)
        (variance 0.0))
    (dotimes (i num-of-samplings)
      (setf lst (cons (monte-carlo-pi num-of-trials) lst)))
    (setf average (/ (reduce #'+ lst) num-of-samplings))
    (setf variance (/ (reduce #'+
                              (mapcar (lambda (x) (expt (- x average) 2)) lst))
                      num-of-samplings))
    (values average variance)))

(defun conf-int-pi (num-of-trials)
  "return the ratio in ~3\sigma."
  (let ((cnter 0)
        (coeff (/ 1.982 (sqrt (1- 100))))
        (average 0.0)
        (variance 0.0))
    (dotimes (i num-of-trials)
      (setf (values average variance) (parameters-pi 100 100))
      (when (and (< (- average (* (sqrt variance) coeff)) pi)
                 (< pi (+ average (* (sqrt variance) coeff))))
        (incf cnter)))
    (/ (* 1.0 cnter) num-of-trials)))
