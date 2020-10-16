module NaiveGradient (naiveGradient) where

import Data.Tuple.Homogenous
import Phylogenetics.Types as Phylo
import qualified Phylogenetics.Likelihood_v2 as V2

eps :: Fractional a => a
eps = 1e-10

naiveGradient
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Topology
  -> NodeMap Double
naiveGradient rate_mx obs bls topo = flip mapWithKey bls $ \node_id bl ->
  let
    bumped_bl = Tuple2 (bl+eps, bl-eps)
    bumped_bls = fmap (\bl' -> insert node_id bl' bls) bumped_bl
    Tuple2 (ll_plus, ll_minus) =
      (V2.logLikelihood rate_mx obs <$> bumped_bls <*> pure topo)
  in
    (ll_plus - ll_minus) / (2 * eps)
