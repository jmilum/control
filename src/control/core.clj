(ns control.core
  (:require
   [apache-commons-matrix.core]
   [clojure.core.matrix :as m]
   [clojure.core.matrix.linear :as ml]
   [control.util :as util]
   [control.state-space :as ss]
   [control.transfer-function :as tf]
   [clojure.math.numeric-tower :as math]))


;(def model-tf {:num [100 20 1] :den [8 2 1] :dt 0})
;(def model-tf (tf/fopdt 1 1 0))
;(def data-step (tf/step model-tf 1 0.1 5))
;(def data-imp (ss/data-step->data-impulse data-step))
;(def model-ss (ss/ho-kalman 32 data-imp))
;(def data-step* (ss/ss->step 50 model-ss))

;(def model-tf (tf/fopdt 1 1 0))
;(def data-step (tf/step model-tf 1 0.01 5))
;(def noise (util/sample-gaussian 0 0.02 500))
;(def data-step-noisy (update data-step :meas #(mapv + noise %)))
;(def data-denoised (ss/filter-svd 2 (:meas data-step-noisy)))
;(def data-imp (ss/data-step->data-impulse {:meas data-denoised :out (into [] (take (count data-denoised) (:out data-step)))}))
;(def model-ss (ss/ho-kalman 1 data-imp))
;(def data-step* (ss/ss->step 500 model-ss))


#_(def model-tf {:num [100 20 1] :den [8 2 1] :dt 0})
#_(def model-tf {:num [1] :den [1 1]})

#_(def model-ss
    (ss/->ss {:A [[0 1] [-8 -2]]
              :B [[0] [1]]
              :C [[92 18]]
              :D [[0]]}))

#_(def resp (lsim-ss model-ss 0.01 (repeat 500 1)))
#_(def resp (lsim-ss model-ss 0.01 (cons 1 (repeat 499 0))))
#_(util/plotxy (:t resp) (:Y resp))

#_(util/plotx (:Y (ss/->step {:model (tf/->ss {:den [1 1 1] :num [5 2 1]}) :tstep 0.1 :tmax 10 :k 1})))
#_(util/plotx (:Y (ss/->impulse {:model (tf/->ss model-tf) :tstep 0.1 :tmax 5 :k 1})))

#_(def G (tf/->ss {:num [3 2 1] :den [-6 5 1]}))

#_(def G (tf/->ss {:num  [1] :den [0 4 5 1]}))
#_(def G (tf/->ss {:num  [1] :den [1 1]}))

#_(util/plotx (:Y (ss/step {:model G :tstep 0.1 :tmax 5 :k 1})))

#_(util/plotx (:Y (ss/lsim-fb G 0.1 (m/matrix [1.473 13.456 18.827]))))
