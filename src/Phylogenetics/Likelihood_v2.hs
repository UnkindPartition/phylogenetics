-- | A different implementation of 'calculatePostOrder'; relies on vector
-- and matrix operations.
module Phylogenetics.Likelihood_v2 where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Storable as VS
import Data.Tuple.Homogenous
import qualified Numeric.LinearAlgebra as Matrix
import Numeric.Log as Log
import Control.Monad.State

import Phylogenetics.Types

-- | The likelihood of a subtree conditional on the state of its root
-- (a.k.a. the «partial likelihood»).
newtype PostOrder = PostOrder { getPostOrder :: Matrix.Vector Double }

-- | The full log-likelihood, summed over all sites
logLikelihood
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Topology
  -> Log Double
logLikelihood rate_mx obs bls tree = product $ do
  site <- [0 .. numOfSites obs - 1]
  return . realToFrac $ likelihood1 rate_mx obs bls site tree

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
    Bin {} -> (VS.sum . getPostOrder) (go topo) / numOfStates rate_mx
  where
  -- Given a (sub)tree topology, return, for each root state,
  -- the probability of observations at the tips
  go :: Topology -> PostOrder
  go = \case
    Leaf leaf_id -> calculateLeafPostOrder rate_mx obs bls site leaf_id
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
calculatePostOrder rate_mx bls subs sub_post_orders =
  let
    Tuple2 (l,r) = do
      sub <- subs
      PostOrder sub_post_order <- sub_post_orders
      let
        bl = bls ! getNodeId sub
        transition_probs = transitionProbabilities rate_mx bl
      pure $ transition_probs Matrix.#> sub_post_order
  in
    PostOrder $ l * r

newtype PreOrder = PreOrder { getPreOrder :: Matrix.Vector Double }

swap :: Tuple2 a -> Tuple2 a
swap (Tuple2 (a,b)) = Tuple2 (b,a)

calculatePreOrder
  :: RateMatrix
  -> BranchLengths
  -> NodeMap PostOrder
  -> PreOrder -- ^ the pre-order of the parent tree
  -> Tuple2 Topology -- ^ the subtrees
  -> Tuple2 PreOrder -- ^ the pre-order traversals of the subtrees
calculatePreOrder rate_mx bls postorders (PreOrder parent_preorder) subs =
  let
    both_transition_probs = do
      sub <- subs
      let bl = bls ! getNodeId sub
      pure $ transitionProbabilities rate_mx bl
  in do
    this_transition_probs <- both_transition_probs
    sibling_transition_probs <- swap both_transition_probs
    sibling <- swap subs
    let
      PostOrder sibling_postorder = postorders ! (getNodeId sibling)
    pure . PreOrder $
      Matrix.tr this_transition_probs Matrix.#>
      (parent_preorder *
        (sibling_transition_probs Matrix.#> sibling_postorder))
