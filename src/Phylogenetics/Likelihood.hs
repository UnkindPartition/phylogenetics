-- | Calculate the likelihood of the tree given the data
module Phylogenetics.Likelihood where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word
import Data.Tuple.Homogenous
import qualified Numeric.LinearAlgebra.Data as Matrix

import Phylogenetics.Types

-- | The likelihood of a subtree conditional on the state of its root.
--
-- If the subtree is a leaf, then it's represented by its character state.
--
-- Otherwise, it's the likelihood conditional on the (unobserved) character
-- state of its root, and is represented as a vector of likelihoods for
-- each state of its root.
data ConditionalLikelihood
  = LeafCL Word8
  | BinCL (VU.Vector Double)

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
  case go topo of
    LeafCL{} -> 1
    BinCL ls -> VU.sum ls / numOfStates rate_mx
  where
  -- Given a (sub)tree topology, return, for each root state,
  -- the probability of observations at the tips
  go :: Topology -> ConditionalLikelihood
  go = \case
    Leaf (NodeId node_id) ->
      LeafCL $ (characters obs IntMap.! node_id) `VU.unsafeIndex` site
    Bin _ sub1 sub2 ->
      let
        subs :: Tuple2 Topology
        subs = tuple2 sub1 sub2
      in
        falg rate_mx bls subs (go <$> subs)

-- | The F-algebra that, given a likelihoods of subtrees calculates the
-- likelihood of the tree
falg
  :: RateMatrix
  -> BranchLengths
  -> Tuple2 Topology
    -- ^ the sub-subtrees (the two immediate children of the current subtree)
  -> Tuple2 ConditionalLikelihood
    -- ^ likelihoods of the sub-subtrees
  -> ConditionalLikelihood
    -- ^ likelihood of the subtree
falg rate_mx bls subs sub_cls = BinCL . VU.fromList $ do
  -- iterate over the possible characters at the root
  root_c <- [0..numOfStates rate_mx - 1]
  return . product $ do
    -- iterate over the two sub-subtrees (applicative do)
    sub <- subs
    sub_cl <- sub_cls
    let
      bl = getBranchLength bls $ getNodeId sub
      transition_probs = transitionProbabilities rate_mx bl
    return $
      case sub_cl of
        LeafCL sub_c -> transition_probs `Matrix.atIndex` (root_c, fromIntegral sub_c)
        BinCL sub_liks -> sum $ do
          -- iterate over the sub-subtree character
          sub_c <- [0..numOfStates rate_mx - 1]
          return $
            VU.unsafeIndex sub_liks (fromIntegral sub_c) *
            (transition_probs `Matrix.atIndex` (root_c, sub_c))
