-- | Calculate the likelihood of the tree given the data
module Phylogenetics.Likelihood where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word
import Data.Tuple.Homogenous
import qualified Numeric.LinearAlgebra.Data as Matrix

import Phylogenetics.Types

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
likelihood1 rate_mx obs bls site = VU.sum . go where
  -- Given a (sub)tree topology, return, for each root state,
  -- the probability of observations at the tips
  go :: Topology -> VU.Vector Double
  go = \case
    Leaf (NodeId node_id) ->
      let
        ch :: Word8 -- this leaf's character at the given site
        ch = (characters obs IntMap.! node_id) `VU.unsafeIndex` site
      in VU.generate (numOfStates rate_mx) (\i -> fromIntegral . fromEnum $ i == fromIntegral ch)
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
  -> Tuple2 (VU.Vector Double)
    -- ^ likelihood of the sub-subtrees (for each of the 4 values of the character)
  -> VU.Vector Double
    -- ^ likelihood of the subtree (for each of the 4 values of the character)
falg rate_mx bls subs sub_liks = VU.fromList $ do
  -- iterate over the possible characters at the root
  root_c <- [0..3]
  return . product $ do
    -- iterate over the two sub-subtrees (applicative do)
    sub <- subs
    sub_lik <- sub_liks
    let
      bl = getBranchLength bls $ getNodeId sub
      transition_probs = transitionProbabilities rate_mx bl
    return . product $ do  
      -- iterate over the sub-subtree character
      sub_c <- [0..3]
      return $
        VU.unsafeIndex sub_lik (fromIntegral sub_c) *
        (transition_probs `Matrix.atIndex` (root_c, sub_c))
