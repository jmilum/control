(ns control.util
  (:require
   [clojure.java.io :as jio]
   [clojure.set :as set]
   [clojure.data.csv :as csv]
   [clojure.string :as str]
   [clojure.math.numeric-tower :as math]
   [clojure.math.combinatorics :as combo]
   [clojure.algo.generic.functor :as f]
   [clojure.core.matrix :as m]
   [com.hypirion.clj-xchart :as chart]
   [clojure.core.matrix.stats :as ms]
   [apache-commons-matrix.core]
   [clojure.core.matrix.linear :as ml]
   [complex.core :as c])
  (:import
   (org.apache.commons.math3.complex Complex)
   (java.io FileNotFoundException BufferedReader)
   (java.util Random)
   (clojure.lang Reflector)
   (org.apache.commons.math3.analysis.solvers LaguerreSolver)))

(set! *warn-on-reflection* true)
(set! *unchecked-math* true)

(m/set-current-implementation :apache-commons)

(defmacro doseq-indexed
  "loops over a set of values, binding index-sym to the 0-based index of each value"
  ([[val-sym values index-sym] & code]
   `(loop [vals# (seq ~values)
           ~index-sym (long 0)]
      (if vals#
        (let [~val-sym (first vals#)]
          ~@code
          (recur (next vals#) (inc ~index-sym)))
        nil))))

(defmacro defmethod*
  "closure over defmethod to bind the dispatch-val with :as"
  [multifn dispatch-val & fn-tail]
  (let [[kw n & body] fn-tail]
    (if (= :as kw)
      `(let [~n ~dispatch-val]
         (defmethod ~multifn ~dispatch-val ~body))
      `(defmethod ~dispatch-val ~fn-tail))))

(defmacro when-let*
  "allow multiple bindings in when-let"
  ([bindings & body]
   (if (seq bindings)
     `(when-let [~(first bindings) ~(second bindings)]
        (when-let* ~(drop 2 bindings) ~@body))
     `(do ~@body))))

(defmacro if-let*
  "allow multiple bindings in if-let"
  ([bindings then]
   `(if-let* ~bindings ~then nil))
  ([bindings then else]
   (if (seq bindings)
     `(if-let [~(first bindings) ~(second bindings)]
        (if-let* ~(drop 2 bindings) ~then ~else)
        ~(if-not (second bindings) else))
     then)))

(defmacro ->map
  "create a map of the values with the names as keywords"
  [& ks]
  (zipmap (map keyword ks) ks))

(defn map-kv [f coll] (reduce-kv (fn [m k v] (assoc m k (f v))) (empty coll) coll))

(defn str->int [s] (when s (Integer/parseInt s)))

(defn numeric? [s]
  (if-let [s (seq s)]
    (let [s (if (= (first s) \-) (next s) s)
          s (drop-while #(Character/isDigit ^Character %) s)
          s (if (= (first s) \.) (next s) s)
          s (drop-while #(Character/isDigit ^Character %) s)]
      (empty? s))))

(defn lazy-file-lines [file]
  (letfn [(helper [^BufferedReader rdr]
            (lazy-seq
             (if-let [line (.readLine ^BufferedReader rdr)]
               (cons line (helper rdr))
               (do (.close rdr) nil))))]
    (helper (jio/reader file))))

(defn load-csv-data
  [file-name separator]
  (with-open [in (jio/reader file-name)]
    (doall (csv/read-csv in :separator separator))))

(defn save-csv-data
  [file-name data]
  (with-open [out (jio/writer file-name)]
    (csv/write-csv out data)))

(defn save-coll
  [file coll]
  (when (seq coll)
    (->> coll
         (interpose \newline)
         (apply str)
         (spit file))))

(defn third [coll] (nth coll 2))
(defn fourth [coll] (nth coll 3))
(defn fifth [coll] (nth coll 4))
(defn sixth [coll] (nth coll 5))
(defn seventh [coll] (nth coll 6))
(defn eighth [coll] (nth coll 7))
(defn ninth [coll] (nth coll 8))
(defn tenth [coll] (nth coll 9))

(defn slurp-from-classpath
  "Slurps a file from the classpath."
  [path]
  (or (some-> path
              jio/resource
              slurp)
      (throw (FileNotFoundException. path))))

(defn get-edn-data
  [file]
  (try
    (read-string (slurp file))
    (catch FileNotFoundException e#
      nil)))

(defn replace-several
  [str & replacements]
  (reduce (fn [s [a b]]
            (str/replace s a b))
          str
          (partition 2 replacements)))

(defn avg-rd [avg i x] (float (+ avg (/ (- x avg) (inc i)))))

(defn reduce-xf [f init]
  (fn [rf]
    (let [acc (volatile! init)]
      (completing (fn [result input] (rf result (vswap! acc f input)))))))

(defn round-dec [n x]
  (Double/parseDouble (format (str "%." n "f") x)))

(defn round-double [n x] (->> x (double) (format (str "%." n "g")) (Double/parseDouble)))

(defn complex? [a] (instance? Complex a))
(defn complex-xs? [xs] (some complex? xs))

(defn inner-prod
  "Inner product of two vectors"
  [xs ys]
  (if (or (complex-xs? xs) (complex-xs? ys))
    (reduce c/+ (map c/* (map c/complex xs)
                     (map c/complex ys)))
    (reduce + (map * xs ys))))

(defn inner-products [xs ys result]
  (if (empty? ys)
    result
    (recur xs (rest ys) (conj result (inner-prod xs ys)))))

(defn convolve
  "Cauchy product of two finite vectors"
  ([xs] xs)
  ([xs ys]
   (if (> (count xs) (count ys))
     (convolve ys xs)                   ;force shorter arg into 1st position
     (let [xs    (reverse xs)
           n     (count xs)
           zeros (repeat (dec n) 0)
           ys    (concat zeros ys)]
       (inner-products xs ys []))))
  ([xs ys & more] (reduce convolve (convolve xs ys) more)))

(defn rand-nums [n] (take n (repeatedly #(rand))))

(def k-min-step 0.02)
(def k-noise-attenuation 3.0)

(defn- norm [vs]
  (Math/sqrt (apply + (map #(Math/pow % 2) vs))))

(defn- clamp [x min max]
  (cond
    (> x max) max
    (< x min) min
    :default x))

(defn low-pass-filter
  ([dt tc adaptive? vals]
   (let [filter-constant (/ dt (+ dt tc))]
     (ref {:filter-constant filter-constant
           :vals            vals
           :adaptive        adaptive?})))
  ([filter input]
   (let [{:keys [filter-constant vals adaptive?]} @filter
         alpha (if adaptive?
                 (let [d (clamp (/ (Math/abs (double (- (norm vals)
                                                        (norm input))))
                                   (- k-min-step 1.0))
                                0.0 1.0)]
                   (+ (/ (* (- 1.0 d) filter-constant)
                         k-noise-attenuation)
                      (* d filter-constant)))
                 filter-constant)
         vals  (map #(+ (* %1 alpha) (* %2 (- 1.0 alpha))) input vals)]
     (:vals (dosync (alter filter assoc :vals vals))))))

(defn high-pass-filter
  ([dt tc adaptive? vals]
   (let [filter-constant (/ tc (+ dt tc))]
     (ref {:filter-constant filter-constant
           :vals            vals
           :raw-vals        vals
           :adaptive        adaptive?})))
  ([filter input]
   (let [{:keys [filter-constant vals raw-vals adaptive?]} @filter
         alpha (if adaptive?
                 (let [d (clamp (/ (Math/abs (double (- (norm vals)
                                                        (norm input))))
                                   (- k-min-step 1.0))
                                0.0 1.0)]
                   (+ (/ (* d filter-constant)
                         k-noise-attenuation)
                      (* (- 1.0 d) filter-constant)))
                 filter-constant)
         vals  (map #(* alpha (- (+ %2 %1) %3)) input vals raw-vals)]
     (:vals (dosync (alter filter assoc :vals vals :raw-vals input))))))

(defn unchunk [s]
  (when (seq s)
    (lazy-seq
     (cons (first s)
           (unchunk (next s))))))

(defn gaussian [mean stdev]
  (let [rng (Random.)]
    (-> (.nextGaussian rng)
        (* stdev)
        (+ mean))))

(defn sample-gaussian [mean stdev n]
  (into [] (repeatedly n #(gaussian mean stdev))))

(defn print-bits [b]
  (let [class-name    (.getName (class b))
        is-byte       (= "java.lang.Byte" class-name)
        num-bits      (Reflector/getStaticField class-name "SIZE")
        format-string (str "~" num-bits "'0b")
        bin-str-fn    #(Reflector/invokeStaticMethod
                        (if is-byte "java.lang.Integer" class-name)
                        "toBinaryString"
                        (to-array [%]))
        bit-string    (if is-byte
                        (str/join (take-last 8 (bin-str-fn (Byte/toUnsignedInt b))))
                        (bin-str-fn b))]
    (str (str/join (repeat (- num-bits (count bit-string)) \0)) bit-string)))

(defn drop-take [coll m]
  (let [start  (get m :start 0)
        amount (get m :amount Integer/MAX_VALUE)
        stop   (min amount (get m :stop Integer/MAX_VALUE))]
    (->> coll (drop start) (take stop))))

(def fib-seq-seq
  ((fn fib [a b]
     (lazy-seq (cons a (fib b (+ a b)))))
   0 1))

(defn pinv [m]
  (let [{:keys [U S V*]}  (ml/svd m)
        Sd (m/diagonal-matrix S)]
    (m/mmul (m/transpose V*) (m/inverse Sd) (m/transpose U))))

(defn msqrt [m]
  (let [{:keys [Q rA]} (ml/eigen m)
        Dsqrt  (m/diagonal-matrix (mapv math/sqrt rA))
        Qinv (m/inverse Q)]
    (m/mmul Q Dsqrt Qinv)))

(defn mpow [m exponent] (apply m/mmul (repeat exponent m)))

(defn hankel [start shift side v]
  (let [v (if (even? (count v)) v (butlast v))
        data (take side (drop shift (partition side 1 (drop start v))))]
    (m/matrix data)))

(defn dehankel [m] (into [] (concat (first (m/columns m)) (rest (last (m/rows m))))))

(defn solve-quadratic
  ([a b c]
   (let [t1 (- 0 b)
         t2 (Math/sqrt (- (* b b) (* 4 a c)))
         t3 (* 2 a)]
     [(/ (- t1 t2) t3)
      (/ (+ t1 t2) t3)])))

(defn symmetric-matrix
  ([data & {:keys [lower] :or {lower true}}]
   (let [n (count data)
         p (int (second (solve-quadratic 1/2 1/2 (- 0 n))))
         mat (m/mutable (m/zero-matrix p p))
         indices (if lower
                   (into [] (for [i (range p) j (range p) :when (<= j i)] [i j]))
                   (into [] (for [i (range p) j (range p) :when (<= i j)] [j i])))]
     (doseq [idx (range n)]
       (let [[i j] (nth indices idx)
             res (nth data idx)]
         (m/mset! mat i j res)
         (m/mset! mat j i res)))
     mat)))

(defn triangular
  ([m & {:keys [upper] :or {upper true}}]
   (let [m (m/mutable m)
         rows (m/row-count m)
         indices (if upper
                   (into [] (for [i (range rows) j (range rows) :when (> i j)] [i j]))
                   (into [] (for [i (range rows) j (range rows) :when (> j i)] [i j])))]
     (doseq [[i j] indices]
       (m/mset! m i j 0))
     m)))

(defn toeplitz [v]
  (symmetric-matrix
   (into []
         (loop [v (rseq v)
                d []]
           (if (nil? v)
             d
             (recur (next v) (concat v d)))))))

(defn peval-1 [v x]
  (reduce-kv
   (fn [sum idx coef] (+ sum (* coef (math/expt x idx))))
   0
   v))

(defn peval
  [m x]
  (/ (peval-1 (:num m) x)
     (peval-1 (:den m) x)))

(defn peval-complex-1 [v x]
  (reduce-kv
   (fn [sum idx coef] (c/+ sum (c/* coef (c/powc x (c/complex idx)))))
   0
   v))

(defn peval-complex
  [m x]
  (let [x' (c/complex x)]
    (c// (peval-complex-1 (:num m) x')
         (peval-complex-1 (:den m) x'))))

(defn pmul-1
  [a b]
  (convolve a b))

(defn pmul
  [a b]
  {:num (convolve (:num a) (:num b))
   :den (convolve (:den a) (:den b))
   :dt (+ (:dt a 0) (:dt b 0))})

(defn tf-pinv [m]
  {:num (:den m)
   :den (:num m)
   :dt 0})

(defn pad-to-n [n x v]
  (into v (repeat (- n (count v)) x)))

(defn padd
  [a b]
  (let [num-term1 (convolve (:num a) (:den b))
        num-term2 (convolve (:den a) (:num b))
        den-term (convolve  (:den a) (:den b))
        max-len (max (count num-term1) (count num-term2))]
    {:num (mapv + (pad-to-n max-len 0 num-term1)
                (pad-to-n max-len 0 num-term2))
     :den (convolve (:den a) (:den b))
     :dt (+ (:dt a 0) (:dt b 0))}))

(defn pdiv [a b]
  (let [b* (tf-pinv b)]
    (pmul a b*)))

(defn psub [a b]
  (let [b* (pmul b {:num [-1] :den [1]})]
    (padd a b*)))

(defn complex-roots [v]
  (-> (LaguerreSolver.)
      (.solveAllComplex (double-array v) 0)
      vec))

(defn combo-prod-sum [v n]
  (reduce
   (fn [acc v]
     (c/+ acc (reduce c/*
                      (c/complex 1)
                      v)))
   (c/complex 0)
   (combo/combinations v n)))

(defn vieta [roots inv-prod order]
  (fn [i]
    (let [sign (if (odd? i) (c/complex -1) (c/complex 1))]
      (c/* sign inv-prod (combo-prod-sum roots (- order i))))))

(defn roots->poly [roots]
  (let [order    (count roots)
        sign     (if (odd? order) (c/complex -1) (c/complex 1))
        inv-prod (c// (reduce c/* roots))
        vieta-fn (vieta roots inv-prod order)
        coef     (mapv c/real-part (map vieta-fn (range 1 order)))]
    (->> (into coef [(c/real-part (c/* sign inv-prod))])
         (into [1]))))

(defn standardize [data]
  (let [mean (ms/mean data)]
    (m/sub data mean)))

(defn normalize [data]
  (let [max (m/maximum data)
        min (m/minimum data)
        apmplitude (+ (math/abs max) (math/abs min))]
    (m/div data apmplitude)))

(defn minima-idx [xs]
  (-> (reduce-kv
       (fn [mins idx [a b c]]
         (when (and (> a b) (< b c))
           (conj mins idx)))
       []
       xs)))

(defn maxima-idx [xs]
  (-> (reduce-kv
       (fn [maxs idx [a b c]]
         (when (and (< a b) (> b c))
           (conj maxs idx)))
       []
       xs)))

(defn min-val-idx [xs]
  (-> (reduce-kv
       (fn [[idx min] k v]
         (if (< min v)
           [idx min]
           [k v]))
       [-1 Double/POSITIVE_INFINITY]
       xs)
      (first)))

(defn max-val-idx [xs]
  (-> (reduce-kv
       (fn [[idx max] k v]
         (if (< v max)
           [idx max]
           [k v]))
       [-1 Double/NEGATIVE_INFINITY]
       xs)
      (first)))

(defn rms [xs]
  (Math/sqrt (/ (reduce + (map #(* % %) xs))
                (count xs))))

(defn plotx [x]
  (chart/view
   (chart/xy-chart
    {"x" {:y x :style {:line-color :red :marker-type :none}}}
    {:title  "Response"
     :x-axis {:title "Time (steps)"}
     :y-axis {:title "% Change"}
     :theme  :matlab})))

(defn plotxy [x y]
  (let [points (min (count x) (count y))
        x'     (take points x)
        y'     (take points y)]
    (chart/view
     (chart/xy-chart
      {"value" {:y y' :x x' :style {:line-color :red :marker-type :none}}}
      {:title  "Response"
       :x-axis {:title "Time (steps)"}
       :y-axis {:title "% Change"}
       :theme  :matlab}))))

(defn plotab [a b]
  (let [points (min (count a) (count b))
        a'     (take points a)
        b'     (take points b)]
    (chart/view
     (chart/xy-chart
      {"a" {:y a' :style {:line-color :red :marker-type :none}}
       "b" {:y b' :style {:line-color :blue :marker-type :none}}}
      {:title  "Response"
       :x-axis {:title "Time (steps)"}
       :y-axis {:title "% Change"}
       :theme  :matlab}))))

(defn pad [n coll val]
  (take n (concat coll (repeat val))))