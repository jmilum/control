(ns control.transfer-function
  (:require
   [clojure.math.numeric-tower :as math]
   [clojure.algo.generic.functor :as f]
   [clojure.core.matrix :as m]
   [clojure.core.reducers :as r]
   [com.hypirion.clj-xchart :as chart]
   [clojure.core.matrix.linear :as ml]
   [clojure.core.matrix.stats :as ms]
   [complex.core :as c]
   [control.util :as util]
   [apache-commons-matrix.core]
   [clojure.core.matrix.protocols :as mp]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* true)

(m/set-current-implementation :apache-commons)

(defn plotm [m k]
  (chart/view
   (chart/xy-chart
    {"y" {:y (k m) :style {:line-color :red :marker-type :none}}}
    {:title  "Response"
     :x-axis {:title "Time (steps)"}
     :y-axis {:title "% Change"}
     :theme  :matlab})))

(defn plot-compare [t y1 y2]
  (let [points (min (count t) (count y1) (count y2))
        t'     (take points t)
        y1'    (take points y1)
        y2'    (take points y2)]
    (chart/view
     (chart/xy-chart
      {"Model" {:x t' :y y2' :style {:marker-color :red :line-style :none}}
       "Data"  {:x t' :y y1' :style {:line-color :blue :marker-type :none}}}
      {:title  "Data vs Model Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab}))))

(defn plot-model [m]
  (let [{:keys [meas t]} m
        points (min (count t) (count meas))
        t      (take points t)
        meas   (take points meas)]
    (chart/view
     (chart/xy-chart
      {"Model" {:x t :y meas :style {:line-color :red :marker-type :none}}}
      {:title  "Model Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab}))))

(defn plot-response [resp]
  (let [{:keys [meas-spt-data out-spt-data meas-dist-data out-dist-data t]} resp
        points (min (count t)
                    (count meas-spt-data)
                    (count out-spt-data)
                    (count meas-dist-data)
                    (count out-dist-data))
        t      (take points t)]
    (chart/view
     (chart/xy-chart
      {"y" {:x     t
            :y     (take points meas-spt-data)
            :style {:line-color :red :marker-type :none}}
       "r" {:x     t
            :y     (repeat points 1)
            :style {:line-color :black :marker-type :none}}}
      {:title  "Setpoint Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab})
     (chart/xy-chart
      {"u" {:x     t
            :y     (take points out-spt-data)
            :style {:line-color :blue :marker-type :none}}}
      {:title  "Setpoint Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab})
     (chart/xy-chart
      {"y" {:x     t
            :y     (take points meas-dist-data)
            :style {:line-color :red :marker-type :none}}
       "r" {:x     t
            :y     (repeat points 0)
            :style {:line-color :black :marker-type :none}}}
      {:title  "Disturbance Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab})
     (chart/xy-chart
      {"u" {:x     t
            :y     (take points out-dist-data)
            :style {:line-color :blue :marker-type :none}}
       "r" {:x     t
            :y     (repeat points 0)
            :style {:line-color :black :marker-type :none}}}
      {:title  "Disturbance Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab}))))

(defn plotr [state]
  (let [{:keys [spt out meas t]} state
        points (count t)
        meas   (if meas (take points meas) (repeat points 0))
        spt    (if spt (take points spt) (repeat points 0))
        out    (if out (take points out) (repeat points 0))]
    (chart/view
     (chart/xy-chart
      {"MEAS" {:x t :y meas :style {:line-color :red :marker-type :none}}
       "SPT"  {:x t :y spt :style {:line-color :blue :marker-type :none}}
       "OUT"  {:x t :y out :style {:line-color :green :marker-type :none}}}
      {:title  "Process Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab}))))

(defn ploty [y ts]
  (let [x (map #(* ts %) (range 0 (count y)))]
    (chart/view
     (chart/xy-chart
      {"y" {:x x :y y :style {:line-color :red :marker-type :none}}}
      {:title  "Response"
       :x-axis {:title "Time (sec)"}
       :y-axis {:title "% Change"}
       :theme  :matlab}))))

(defn huddleston
  [fs t]
  (let [sigma (+ 0.01 (apply max (map c/real-part (util/complex-roots (:den fs)))))
        nint  500
        delta 0.01
        k     (/ (* (Math/exp (* sigma t)) delta 0.5) 3.14159)
        f     (fn [i]
                  (c/real-part
                   (c/* (c/exp (c/complex 0 (* t i)))
                        (util/peval-complex fs (c/complex sigma i)))))
        rf    (fn
               ([] 0)
               ([sum i] (+ sum (f i) (f (+ i delta)))))]
    (* k (r/fold + rf (vec (range 0 nint delta))))))

(defn zakian [fs t]
  (let [k [[(c/complex 12.83767675 1.666063445) (c/complex -36902.08210 196990.4257)]
           [(c/complex 12.22613209 5.012718792) (c/complex 61277.02524 -95408.62551)]
           [(c/complex 10.93430308 8.409673116) (c/complex -28916.56288 18169.18531)]
           [(c/complex 8.776434715 11.92185389) (c/complex 4655.361138 -1.901528642)]
           [(c/complex 5.225453361 15.72952905) (c/complex -118.7414011 -141.3036911)]]]
    (-> (reduce
         (fn [acc [k1 k2]]
             (+ acc
                (c/real-part
                 (c/* k2
                      (util/peval-complex fs (c// k1 t))))))
         0
         k)
        (* 2)
        (/ t))))

(defn least-squares
  [data dt order]
  (let [{:keys [out meas t tstep]} data
        shift     (/ dt tstep)
        y-shifted (drop shift meas)
        a-data    (partition order 1 (concat (repeat order 0) (map - y-shifted)))
        b-data    (partition order 1 (concat (repeat (dec order) 0) out))
        y-data    (partition 1 1 y-shifted)
        rows      (min (count a-data) (count b-data) (count y-data))
        A         (m/clone (m/join-along 1 (m/matrix (take rows a-data)) (m/matrix (take rows b-data))))
        Y         (m/to-vector (take rows y-data))
        X         (ml/least-squares A Y)
        params    (m/to-vector X)
        num       (mapv #(util/round-double 6 %) (vec (drop order params)))
        den       (mapv #(util/round-double 6 %) (conj (vec (take order params)) 1))
        Y*        (m/mmul A X)
        SSE       (m/mget (ms/sum-of-squares (m/sub Y Y*)))
        AIC       (+ (* 2 (inc (+ order order))) (* rows (Math/log SSE)))]
    (util/->map AIC num den tstep dt)))

(defn model-fitness
  [data dt order]
  (:AIC (least-squares data dt order)))

(defn fit-dt [data order start end]
  (->> (range start end (:tstep data))
       (map #(model-fitness data % order))))

(defn z-root->s-root
  [tstep root]
  (let [ln-root (c/log root)]
    (if (= ln-root (c/complex 0 0))
      ln-root
      (c// tstep (c/log root)))))

(defn s-root->z-root [tstep root] (c/exp (c/* root tstep)))

(defn complex-roots->poly [roots]
  (->> roots
       (map c/-)
       (map (fn [x]
                (if (= x (c/complex 0 0))
                  [x (c/complex 1 0)]
                  [(c/complex 1 0) x])))
       (apply util/convolve)
       (mapv #(c/real-part %))))

(defn round-complex [a]
  (let [r (c/real-part a)
        i (c/imaginary-part a)]
    (c/complex (util/round-double 4 r) (util/round-double 4 i))))

(defn z->s [m]
  (let [{:keys [num den dt tstep]} m
        num*     (mapv #(util/round-double 6 %) num)
        den*     (mapv #(util/round-double 6 %) den)
        den-eval (util/peval-1 den* 1)
        num-eval (util/peval-1 num* 1)
        k        (if (= 0.0 (util/round-double 6 den-eval))
                   (/ (util/peval-1 num 1) (math/expt tstep (:order m)))
                   (util/peval m 1))
        s-den    (->> (util/complex-roots den*)
                      (map round-complex)
                      (mapv #(z-root->s-root tstep %))
                      (complex-roots->poly)
                      (mapv #(util/round-double 4 %)))]
    {:dt  dt
     :den s-den
     :num [k]
     :k   k}))

(defn s->z [m]
  (let [{:keys [num den dt tstep]} m
        num*     (mapv #(util/round-double 6 %) num)
        den*     (mapv #(util/round-double 6 %) den)
        den-eval (util/peval-1 den* 1)
        num-eval (util/peval-1 num* 1)
        k        (if (= 0.0 (util/round-double 6 den-eval))
                   (/ (util/peval-1 num 1) (math/expt tstep (:order m)))
                   (util/peval m 1))
        s-den    (->> (util/complex-roots den*)
                      (map round-complex)
                      (mapv #(s-root->z-root tstep %))
                      (complex-roots->poly)
                      (mapv #(util/round-double 4 %)))]
    {:dt  dt
     :den s-den
     :num [k]
     :k   k}))

(defn data->z [data dt order]
  (let [tmax    (last (:t data))
        tsamp   (- (second (:t data)) (first (:t data)))
        z-model (least-squares data dt order)]
    (merge z-model {:order order})))

(defn step
  [model k tstep tmax]
  (let [{:keys [num den dt]} (merge {:dt 0} model)
        dt*    (if (< dt 0) 0 dt)
        f      (util/pmul model {:num [k] :den [0 1]})
        t      (mapv #(util/round-double 5 %) (butlast (range 0 tmax tstep)))
        shift  (/ dt* tstep)
        ty-end (dec (- (/ tmax tstep) shift))
        ty     (subvec t 1 ty-end)
        meas   (into (vec (repeat shift 0)) (mapv #(zakian f %) ty))
        out    (into [] (repeat (dec (count meas)) k))]
    (merge model (util/->map meas out t tmax tstep))))

(defn pid [p i d]
  (let [Kc  (/ 100 p)
        Ti  (* i 60)
        Td  (* d 60)
        n   10
        tau (/ Td n)]
    {:num [Kc (* Kc (+ Ti (/ Td n))) (* Kc (+ (* Ti Td) (* Ti tau)))]
     :den [0 Ti (* Ti tau)]}))

(defn add-noise [data k mean stdev]
  (assoc data k (m/add (k data) (util/sample-gaussian mean stdev (m/row-count (k data))))))

(defn fopdt [k tau dt] {:num [k] :den [1 tau] :dt dt :k k :tau tau})
(defn sopdt [k tau zeta dt] {:num [k] :den [1 (* 2 tau zeta) (* tau tau)] :dt dt})
(defn intdt [k dt] {:num [k] :den [0 1] :k k :dt dt})
(defn dintdt [k dt] {:num [k] :den [0 1 1] :dt dt})

(defn pband [model tc tau] (/ (* 100 (:k model) (+ tc (:dt model))) tau))

(defn get-tuning-2+ [model tc]
  (let [roots (util/complex-roots (:den model))
        taus  (map #(c// -1.0 %) roots)
        tau+  (c/real-part (reduce c/+ taus))
        tau*  (c/real-part (reduce c/* taus))
        PBAND (pband model tc tau+)
        INT   (/ tau+ 60)
        DERIV (/ tau* tau+ 60)]
    (util/->map PBAND INT DERIV)))

(defn get-tuning-1 [model tc]
  (let [tau   (second (:den model))
        PBAND (pband model tc tau)
        INT   (/ tau 60)
        DERIV 0]
    (util/->map PBAND INT DERIV)))

(defn tune [model tc]
  (if
   (< (count (:den model)) 2)
    (get-tuning-1 model tc)
    (get-tuning-2+ model tc)))

(defn tune-int [model]
  (let [dt    (:dt model)
        k     (:k model)
        PBAND (* 200 k dt)
        INT   (* 6.8 dt)
        DERIV 0]
    (util/->map PBAND INT DERIV)))

(defn process-response
  [process controller step-size tstep tmax]
  (let [open-loop      (util/pmul controller process)
        den            (util/padd open-loop {:num [1] :den [1] :dt 0})
        meas-spt       (util/pdiv open-loop den)
        out-spt        (util/pdiv controller den)
        meas-dist      (assoc (util/pdiv process den) :dt 0)
        out-dist       (assoc (util/pmul meas-spt {:num [-1] :den [1] :dt 0}) :dt 0)
        meas-spt-data  (:y (step meas-spt step-size tstep tmax))
        out-spt-data   (:y (step out-spt step-size tstep tmax))
        meas-dist-data (:y (step meas-dist step-size tstep tmax))
        out-dist-data  (:y (step out-dist step-size tstep tmax))
        t              (range 0 tmax tstep)]
    (util/->map open-loop den meas-spt out-spt meas-dist out-dist meas-spt-data
                out-spt-data meas-dist-data out-dist-data t)))

(def k 0.000195)
(def tau 60)
(def dt 300)
(def zeta 1.06)
(def tstep 0.1)
(def tmax 60)
(def step-size 1)

;(def model (fopdt k tau dt))
;(def model (sopdt k tau zeta dt))
;(def model (intdt k dt))
(def model {:num [5 5 1] :den [2 6.5 5 1.65 1] :dt 0})

#_(def data (step model step-size tstep tmax))
#_(def data (add-noise (step model step-size tsamp tmax) :y 0 0.0))
#_(plot-model data)

;(def order 1)
#_(plotx (fit-dt data order))
(def dt-est 5)

;(def z-model (data->z data dt-est order))
;(def s-model (z->s z-model))
;(def s-data (step s-model step-size tsamp tmax))
#_(plot-compare (:t data) (:y data) (:y s-data))

;(def agro 1)
;(def use-deriv 0)
#_(def tuning (tune s-model 1))
;(def tuning (tune model 1))
;(def tuning (tune-int model))

;(def resp (process-response
;            model
;            ;(pid (/ (:PBAND tuning) agro) (:INT tuning) (* use-deriv (:DERIV tuning)))
;            (pid 4 30 0)
;            9
;            tsamp
;            tmax))

(defn slice-response [state options]
  (let [{:keys [spt out meas t]} state
        options (f/fmap #(/ % (:tstep state)) options)
        slice   #(util/drop-take % options)]
    (-> state
        (update :meas slice)
        (update :spt slice)
        (update :out slice)
        (update :t slice))))

(def nth* (fnil nth '()))

(defn fop [state]
  (let [{:keys [k tstep tau dt out meas bias-out]} state
        meas0 (nth* meas 0 0)
        dt    (if dt (/ dt tstep) 0)
        out0  (nth* out dt 0)
        out   (if bias-out (+ out0 bias-out) out0)
        a     (- 1 (/ tstep tau))
        b     (/ (* k tstep) tau)
        val   (max 0 (min 100 (+ (* a meas0) (* b out))))]
    (update state :meas #(conj % val))))

(defn sop [state]
  (let [{:keys [k zeta tstep tau out meas dt bias-out]} state
        dt    (if dt (/ dt tstep) 0)
        meas0 (nth* meas 0 0)
        meas1 (nth* meas 1 0)
        out0  (nth* out dt 0)
        out   (if bias-out (+ out0 bias-out) out0)
        a     (/ (* k tstep tstep) (* tau tau))
        b     (/ (* 2 (- tau (* zeta tstep))) tau)
        c     (/ (- (* 2 zeta tau tstep) (* tau tau) (* tstep tstep)) (* tau tau))
        val   (->> (+ (* a out)
                      (* b meas0)
                      (* c meas1))
                   (min 100)
                   (max 0))]
    (update state :meas #(conj % val))))

(defn integrating [state]
  (let [{:keys [k out meas dt tstep]} state
        dt    (/ dt tstep)
        meas0 (nth* meas 0 0)
        out0  (nth* out dt 0)
        val   (->> (+ meas0 (* k out0))
                   (min 100)
                   (max 0))]
    (assoc state :meas val)))

(defn butterworth [state]
  (let [{:keys [kd int deriv tstep meas measb]} state
        meas0  (nth* meas 0 0)
        measb0 (nth* measb 0 0)
        measb1 (nth* measb 1 0)
        kd     (if kd kd 10)
        tau    (/ (* kd (+ (/ int) (/ deriv))))
        int    (/ (* 60 int) tstep)
        deriv  (/ (* 60 deriv) tstep)
        a      (/ (* tstep tstep) (* 0.5 tau tau))
        b      (/ (* 2 (- (* 0.5 tau tau) (* tau tstep))) (* tau tau))
        c      (/ (- (* tau tstep) (* 0.5 tau tau) (* tstep tstep)) (* 0.5 tau tau))
        val    (->> (+ (* a meas0)
                       (* b measb0)
                       (* c measb1))
                    (min 100)
                    (max 0))]
    (update state :measb #(conj % val))))

(defn high-pass [state]
  (let [{:keys [tstep meas measb int deriv kd] :or {kd 10}} state
        int    (* 60 int)
        deriv  (if deriv (* 60 deriv))
        tau    (/ (* int deriv) (* kd (+ deriv int)))
        meas0  (nth* meas 0 0)
        meas1  (nth* meas 1 0)
        measb0 (nth* measb 0 0)
        a      (if (zero? deriv) 0 (- 1 (/ tstep tau)))
        val    (+ meas0 (* -1 meas1) (* a measb0))]
    (update state :measb #(conj % val))))

(defn low-pass [state]
  (let [{:keys [tlpass tstep meas measb]} state
        meas0  (nth* meas 0 0)
        measb0 (nth* meas 0 0)
        a      (- 1 (/ tstep tlpass))
        b      (/ tstep tlpass)
        val    (max 0 (min 100 (+ (* a measb0) (* b meas0))))]
    (update state :measb #(conj % val))))

(defn pid [state]
  (let [{:keys [pband spllag int deriv spt out meas measb tstep]} state
        meas0  (nth* meas 0 0)
        meas1  (nth* meas 1 0)
        meas2  (nth* meas 2 0)
        measb0 (nth* measb 0 0)
        measb1 (nth* measb 1 0)
        measb2 (nth* measb 2 0)
        out0   (nth* out 0 0)
        spt0   (nth* spt 0 0)
        spt1   (nth* spt 1 0)
        int    (* 60 int)
        deriv  (if deriv (* 60 deriv) 0)
        k      (/ 100 pband)
        a      (/ (* tstep k) int)
        b      (/ (* -1 k deriv) int)
        c      (/ (* -1 k deriv) tstep)
        e0*    (- (* spt0 spllag) meas0)
        e1     (- spt1 meas1)
        e1*    (- (* spt1 spllag) meas1)
        c1     (- measb0 measb1)
        c2     (+ measb0 (* -2 measb1) measb2)
        val    (max 0 (min 100 (+ out0 (* a e1) (* k e0*) (* k -1 e1*) (* b c1) (* c c2))))]
    (-> state
        high-pass
        (update :out #(conj % val)))))

(defn pi [state]
  (let [{:keys [pband spllag int spt out meas tstep]} state
        meas0 (nth* meas 0 0)
        meas1 (nth* meas 1 0)
        out0  (nth* out 0 0)
        spt0  (nth* spt 0 0)
        spt1  (nth* spt 1 0)
        int   (* 60 int)
        k     (/ 100 pband)
        a     (/ (* tstep k) int)
        e0*   (- (* spt0 spllag) meas0)
        e1    (- spt1 meas1)
        e1*   (- (* spt1 spllag) meas1)
        val   (max 0 (min 100 (+ out0 (* a e1) (* k e0*) (* k -1 e1*))))]
    (update state :out #(conj % val))))

(defn time-inc [state]
  (let [{:keys [tstep t]} state
        t0 (nth* t 0 0)
        tn (+ t0 tstep)]
    (update state :t #(conj % tn))))

(defn simulate [state f]
  (loop [state state]
    (let [{:keys [t tmax dynamics]} state
          t0      (nth* t 0 0)
          process (comp dynamics f time-inc)]
      (if (<= t0 tmax)
        (recur (process state))
        (-> state
            (update :meas (comp vec reverse))
            (update :spt (comp vec reverse))
            (update :out (comp vec reverse))
            (update :t (comp vec reverse)))))))

(defn impulse-out [k]
  (fn [state]
      (if (nil? (:out state))
        (assoc state :out k)
        (assoc state :out 0))))

(defn assoc-param [param f]
  (fn [state]
      (let [coll (get state param)]
        (assoc state param (conj (pop coll) (f (nth* coll 0 0)))))))

(defn update-param [param f]
  (fn [state]
      (let [coll (get state param)]
        (update state param (conj (pop coll) (f (nth* coll 0 0)))))))

(defn conj-const-param [param val]
  (fn [state]
      (update state param #(conj % val))))

(defn prbs-gen [param delta-t low high]
  (fn [state]
      (let [{:keys [tstep t]} state
            val (cond
                  (<= tstep (mod (nth* t 0 0) delta-t)) (nth* (get state param) 0 0)
                  (= 0 (rand-int 2)) low
                  :else high)]
        (update state param #(conj % val)))))

(defn s->2nd-order [m]
  (let [{:keys [k den dt]} m
        tau  (math/sqrt (nth den 2))
        zeta (/ (nth den 1) (* 2 tau))]
    (util/->map dt k tau zeta)))

(defn s->1st-order [m]
  (let [{:keys [k den dt]} m
        tau (nth den 1)]
    (util/->map dt k tau)))

(defn error-itae [state]
  (-> (m/to-vector (:spt state))
      (m/sub (m/to-vector (:meas state)))
      (mp/abs)
      (m/mul (m/to-vector (:t state)))
      (ms/sum)))

(defn error-iae [state]
  (-> (m/to-vector (:spt state))
      (m/sub (m/to-vector (:meas state)))
      (mp/abs)
      (ms/sum)))

(defn error-ise [state]
  (-> (m/to-vector (:spt state))
      (m/sub (m/to-vector (:meas state)))
      (ms/sum-of-squares)))

(defn error-itse [state]
  (-> (m/to-vector (:spt state))
      (m/sub (m/to-vector (:meas state)))
      (m/square)
      (m/mul (m/to-vector (:t state)))
      (ms/sum)))

(defn error-mse [state]
  (-> (m/to-vector (:spt state))
      (m/sub (m/to-vector (:meas state)))
      (ms/sum-of-squares)
      (/ (count (:spt state)))))

(defn loop-exercise [m]
  (let [{:keys [spt bias-out switch-time] :or {spt 50 bias-out 0 switch-time 0}} m]
    (fn [state]
        (let [t    (nth* (:t state) 0 0)
              bias (if (<= t switch-time) 0 bias-out)]
          (-> state
              (assoc :bias-out bias)
              (update :spt #(conj % spt)))))))

(defn parallel->serial [m]
  (let [{:keys [kc tau-i tau-d] :or {tau-d 0}} m
        fc    (Math/sqrt (- 1 (* 4 (/ tau-d tau-i))))
        kc    (* 0.5 kc (inc fc))
        tau-i (* 0.5 tau-i (inc fc))
        tau-d (* 0.5 kc (- 1 fc))]
    (util/->map fc kc tau-i tau-d)))

(defn cohen-coon-pi [m]
  (let [{:keys [k tau dt]} m
        kc    (* (/ k) (/ tau dt) (+ 0.9 (/ dt (* 12 tau))))
        tau-i (* dt (/ (+ 30 (/ (* 3 dt) tau)) (+ 9 (/ (* 20 dt) tau))))]
    (util/->map kc tau-i)))

(defn itae-spt-pid [m]
  (let [{:keys [k tau dt]} m
        kc    (* (/ 0.965 k) (Math/pow (/ dt tau) -0.855))
        tau-i (/ tau (+ 0.796 (* -0.147 (/ dt tau))))
        tau-d (* 0.308 tau (Math/pow (/ dt tau) 0.9292))]
    (util/->map kc tau-i tau-d)))

(defn chr-0%-spt-pi [m]
  (let [{:keys [k tau dt]} m
        kc    (/ (* 0.35 tau) (* k dt))
        tau-i (* 1.2 tau)]
    (util/->map kc tau-i)))

(defn chr-20%-spt-pi [m]
  (let [{:keys [k tau dt]} m
        kc    (/ (* 0.6 tau) (* k dt))
        tau-i tau]
    (util/->map kc tau-i)))

(defn s->model [m]
  {:k   (:k m)
   :tau (last (:den m))
   :dt  (:dt m)})

(defn serial->pida [m]
  (let [{:keys [kc tau-i tau-d] :or {tau-d 0}} m]
    {:pband (/ 100 kc)
     :int   (/ tau-i 60)
     :deriv (/ tau-d 60)}))

(defn step-info [m]
  (let [{:keys [t meas]} m
        peak-index      (util/max-val-idx meas)
        peak-value      (nth meas peak-index)
        peak-time       (nth t peak-index)
        meas-head       (util/drop-take meas {:stop peak-index})
        meas-tail       (util/drop-take meas {:start peak-index})
        settling-min    (apply min meas-tail)
        t-tail          (util/drop-take t {:start peak-index})
        meas-final      (ms/mean (util/drop-take meas-tail {:start (Math/floor (/ (count meas-tail) 2))}))
        rise-low-index  (util/min-val-idx (m/to-vector (mp/abs (m/sub meas-head (* 0.1 meas-final)))))
        rise-high-index (util/min-val-idx (m/to-vector (mp/abs (m/sub meas-head (* 0.9 meas-final)))))
        rise-time       (- (nth t rise-high-index) (nth t rise-low-index))
        settling-index  (util/min-val-idx (mp/abs (m/sub (m/sub meas-tail meas-final) (* 0.02 peak-value))))
        settling-time   (nth t-tail settling-index)
        overshoot       (* 100 (/ (- peak-value meas-final) meas-final))
        undershoot      (* 100 (/ (- meas-final settling-min) meas-final))]
    (util/->map peak-value peak-time rise-time settling-time overshoot undershoot)))

(defn fitness [m weights]
  (let [{:keys [rise-time settling-time overshoot undershoot]} (step-info m)
        itae (error-itae m)
        {:keys [a b c d e]} weights]
    (+ (* a itae) (* b overshoot) (* c undershoot) (* d rise-time) (* e settling-time))))

#_(def process {:pband    200 :int 0.2 :deriv 0.25 :spllag 1 :k 1.25 :tau 10 :zeta 0.9 :dt 0
                :dynamics sop :tstep 0.1 :tmax 150 :kd 10})
#_(def process-response (loop-gen process (comp pid (conj-const-param :spt 50))))

(def process {:pband    105 :int 0.75 :deriv 0.1 :spllag 1 :k 1 :tau 25 :zeta 0.9 :dt 10
              :dynamics sop :tstep 0.1 :tmax 240})

#_(def process-response (simulate process (comp pid (loop-exercise {:out 20 :bias-out 0 :switch-time 0}))))

#_(fitness process-response {:a 0.5 :b 10 :c 5 :d 2 :e 10})

#_(plotr process-response)
#_(error-ise process-response)

#_(step-info process-response)
#_(error-itae process-response)

#_(plotx (fit-dt (resp->data process-response) 2 0 15))
#_(-> process-response
      (slice-response {:stop 240})
      resp->data
      (data->z 10 1)
      z->s
      s->model
      ;cohen-coon-pi
      ;itae-spt-pid
      ;chr-0%-spt-pi
      ;chr-20%-spt-pi
      parallel->serial
      serial->pida)

#_(step-info process-response)

(defn construct-A [m]
  (let [n     (dec (count (:den m)))
        ones  (m/identity-matrix (dec n))
        zeros (m/column-matrix (repeat (dec n) 0))
        data  (m/row-matrix (map unchecked-negate (butlast (:den m))))]
    (if (< n 2)
      data
      (m/join-along 0 (m/join-along 1 zeros ones) data))))

(defn construct-B [m]
  (let [n (dec (count (:den m)))]
    (m/column-matrix (reverse (util/pad n [1] 0)))))

;(defn construct-C [m]
;  (let [{:keys [num den]} m
;        n (dec (count den))]
;    (m/row-matrix (util/pad n num 0))))

(defn construct-C [m]
  (let [{:keys [num den]} m
        n  (count den)
        b' (util/pad n num 0)
        b0 (last b')
        b  (butlast b')
        a  (butlast den)]
    (m/row-matrix (mapv #(- %1 (* %2 b0)) b a))))

(defn construct-D [m]
  (let [{:keys [num den]} m
        n (count den)]
    (m/matrix [[(last (util/pad n num 0))]])))

(defn ->ss [m]
  {:A (construct-A m)
   :B (construct-B m)
   :C (construct-C m)
   :D (construct-D m)})
