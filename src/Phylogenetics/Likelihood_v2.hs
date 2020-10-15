-- | A different implementation of 'calculatePostOrder'; relies on vector
-- and matrix operations.
module Phylogenetics.Likelihood_v2 where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Storable as VS
import Data.Tuple.Homogenous
import qualified Numeric.LinearAlgebra as Matrix

import Phylogenetics.Types

-- | The likelihood of a subtree conditional on the state of its root
-- (a.k.a. the «partial likelihood»).
--
-- If the subtree is a leaf, then it's represented by its character state.
--
-- Otherwise, it's the likelihood conditional on the (unobserved) character
-- state of its root, and is represented as a vector of likelihoods for
-- each state of its root.
type PostOrder = Matrix.Vector Double

-- | The full log-likelihood, summed over all sites
logLikelihood
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Topology
  -> Double
logLikelihood rate_mx obs bls tree = sum $ do
  site <- [0 .. numOfSites obs - 1]
  return . log $ likelihood1 rate_mx obs bls site tree

-- | Likelihood for a single site
likelihood1
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Int -- ^ the index of the site
  -> Topology
  -> Double
likelihood1 rate_mx obs bls site topo =
  case topo of
    Leaf {} -> 1
    Bin {} -> VS.sum (go topo) / numOfStates rate_mx
  where
  -- Given a (sub)tree topology, return, for each root state,
  -- the probability of observations at the tips
  go :: Topology -> PostOrder
  go = \case
    Leaf (NodeId node_id) ->
      let ch = (characters obs IntMap.! node_id) `VU.unsafeIndex` site
      in VS.generate (numOfStates rate_mx) $ \i -> if i == fromIntegral ch then 1 else 0
    Bin _ sub1 sub2 ->
      let
        subs :: Tuple2 Topology
        subs = tuple2 sub1 sub2
      in
        calculatePostOrder rate_mx bls subs (go <$> subs)

calculatePostOrder
  :: RateMatrix
  -> BranchLengths
  -> Tuple2 Topology
    -- ^ the subtrees
  -> Tuple2 PostOrder
    -- ^ likelihoods of the subtrees
  -> PostOrder
    -- ^ likelihood of the subtree
calculatePostOrder rate_mx bls subs sub_cls =
  let
    Tuple2 (l,r) = do
      sub <- subs
      sub_cl <- sub_cls
      let
        bl = getBranchLength bls $ getNodeId sub
        transition_probs = transitionProbabilities rate_mx bl
      pure $ transition_probs Matrix.#> sub_cl
  in
    l * r
