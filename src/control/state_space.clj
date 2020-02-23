(ns control.state-space
  (:require
   [clojure.math.numeric-tower :as math]
   [clojure.algo.generic.functor :as f]
   [clojure.core.reducers :as r]
   [com.hypirion.clj-xchart :as chart]
   [complex.core :as c]
   [apache-commons-matrix.core]
   [control.util :as util]
   [clojure.core.matrix :as m]
   [clojure.core.matrix.linear :as ml]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* true)

(m/set-current-implementation :apache-commons)

(defn filter-svd [n data]
  (let [H (util/hankel 0 0 (int (math/floor (* 0.5 (count data)))) data)
        {:keys [U S V*]} (ml/svd H)
        rows (m/row-count S)
        Sd (m/diagonal-matrix (concat (take n S) (take (- rows n) (repeat 0))))
        H' (m/mmul U Sd V*)]
    (util/dehankel H')))

(defn data-step->data-impulse [data]
  (let [{:keys [meas out]} data]
    (->> out
         (util/toeplitz)
         (util/triangular)
         (util/pinv)
         (m/mmul (take (count out) meas))
         (into []))))

;; impulse response from markov parameters
;(defn ->impulse
;  ([n model]
;   (->impulse n model 0))
;  ([n model i]
;   (let [{:keys [A B C D]} model
;         term (cond
;                (== i 0) (m/mget D 0 0)
;                (== i 1) (m/mget (m/mmul C B) 0 0)
;                (> i 1) (m/mget (m/mmul C (util/mpow A (dec i)) B) 0 0))]
;     (take n (lazy-seq (cons term (->impulse n model (inc i))))))))

(defn ho-kalman [n data]
  (let
   [H0 (util/hankel 1 0 4 data)
    H1 (util/hankel 1 1 4 data)
    {:keys [U S V*]} (ml/svd H0)
    V (m/transpose V*)
    Us (m/transpose (m/matrix (take n (m/columns U))))
    Ss (m/diagonal-matrix (take n S))
    Ssq  (util/msqrt Ss)
    Vs (m/matrix (take n (m/columns V)))
    Ok (m/mmul Us Ssq)
    Cl (m/mmul Ssq Vs)
    Ahat (m/to-nested-vectors (m/mmul (util/pinv Ok) H1 (util/pinv Cl)))
    Bhat (m/to-nested-vectors (m/transpose (m/matrix (take 1 (m/columns Cl)))))
    Chat (m/to-nested-vectors (take 1 (m/rows Ok)))
    Dhat (m/to-nested-vectors (first data))]
    {:A Ahat :B Bhat :C Chat :D Dhat}))

(defn ->ss [m]
  (-> m
      (update :A m/matrix)
      (update :B m/matrix)
      (update :C m/matrix)
      (update :D m/matrix)))

(defn continous->discrete [m tstep]
  (let [{:keys [A B C D]} m
        a (/ 2 tstep)
        I (m/identity-matrix (m/row-count A))
        aI (m/mmul a I)
        aIA (util/pinv (m/sub aI A))
        Ad (m/mmul aIA (m/add aI A))
        Bd (m/mmul aIA B)
        Cd (m/mmul C (m/add Ad I))
        Dd (m/add (m/mmul C Bd) D)]
    (merge m (util/->map Ad Bd Cd Dd tstep))))

(defn ssd-eval [model Un]
  (let [{:keys [Ad Bd Cd Dd X t tstep]} model
        Xn (last X)
        tn (last t)
        Yn (m/add (m/mmul Cd Xn) (m/mmul Dd Un))
        Xn+1 (m/add (m/mmul Ad Xn) (m/mmul Bd Un))]
    (-> model
        (update :X conj Xn+1)
        (update :Y conj Yn)
        (update :t conj (+ tn tstep)))))

(defn lsim [model tstep U]
  (let [order (m/row-count (:A model))
        m (merge model {:t [0] :X [(m/matrix (vec (repeat order [0])))]})
        md (continous->discrete m tstep)]
    (-> (reduce ssd-eval md U)
        (update :Y (partial map #(m/mget % 0 0)))
        (update :Y reverse)
        (update :Y vec)
        (dissoc :X))))

(defn lsim->fn [f]
  (fn [{:keys [model tstep tmax k]}]
    (let [n (/ tmax tstep)
          U (f n k)]
      (lsim model tstep U))))

(def step (lsim->fn #(repeat %1 %2)))
(def impulse (lsim->fn #(cons %2 (repeat (dec %1) 0))))

(defn ssd-eval-fb [model _]
  (let [{:keys [Ad Bd Cd X t tstep]} model
        Xn (last X)
        tn (last t)
        Yn (m/mmul Cd Xn)
        Xn+1 (m/sub (m/mmul Ad Xn) Bd)]
    (-> model
        (update :X conj Xn+1)
        (update :Y conj Yn)
        (update :t conj (+ tn tstep)))))

(defn lsim-fb [model tstep K]
  (let [{:keys [A B]} model
        Acl (m/sub A (m/mmul K B))
        order (m/row-count (:A model))
        m (merge model {:A Acl :t [0] :X [(m/matrix (vec (repeat order [0])))]})
        md (continous->discrete m tstep)]
    (-> (reduce ssd-eval-fb md (range 50))
        (update :Y (partial map #(m/mget % 0 0)))
        (update :Y reverse)
        (update :Y vec)
        (dissoc :X))))
