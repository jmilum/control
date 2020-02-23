(ns user
  (:require
    [clojure.set :as set]
    [clojure.string :as str]
    [clojure.math.numeric-tower :as math]
    [clojure.math.combinatorics :as combo]
    [clojure.algo.generic.functor :as f]
    [com.hypirion.clj-xchart :as chart]
    [clojure.pprint :refer [pprint pp]]
    [clojure.repl :refer [apropos dir doc find-doc pst source]]
    [clojure.tools.namespace.repl :refer [refresh refresh-all]]
    [criterium.core :refer [bench quick-bench]]
    [clojure.edn :as edn]
    [control.core :refer :all]
    [control.util :as util]
    [apache-commons-matrix.core]
    [clojure.core.matrix :as m]
    [clojure.core.matrix.linear :as ml]
    [control.state-space :as ss]
    [control.transfer-function :as tf]
    [clojure.math.numeric-tower :as math]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* true)

(m/set-current-implementation :apache-commons)
