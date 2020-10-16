-- | A different implementation of 'calculatePostOrder'; relies on vector
-- and matrix operations.
module Phylogenetics.Likelihood_v2 where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Storable as VS
import Data.Tuple.Homogenous
import qualified Numeric.LinearAlgebra as Matrix
import Control.Monad.State
import Data.List (foldl')

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
    Bin {} -> (VS.sum . getPostOrder) (go topo) / numOfStates rate_mx
  where
  -- Given a (sub)tree topology, return, for each root state,
  -- the probability of observations at the tips
  go :: Topology -> PostOrder
  go = \case
    Leaf leaf_id -> calculateLeafPostOrder rate_mx obs site leaf_id
    Bin _ sub1 sub2 ->
      let
        subs :: Tuple2 Topology
        subs = tuple2 sub1 sub2
      in
        calculatePostOrder rate_mx bls subs (go <$> subs)

-- | Calculate the log-likelihood and the gradient of the log-likelihood
-- w.r.t. all branch lengths
gradient
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Topology
  -> (Double, NodeMap Double) -- ^ log-likelihood and its gradient
gradient rate_mx obs bls tree =
  foldl' (\(l,g) (l1,g1) -> ((,) $! (l+l1)) $! (g+g1)) (0, mempty) $ do
    site <- [0 .. numOfSites obs - 1]
    return $ gradient1 rate_mx obs bls site tree

-- | Log-likelihood and gradient for a single site
gradient1
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Int -- ^ the index of the site
  -> Topology
  -> (Double, NodeMap Double)
gradient1 rate_mx obs bls site topo =
  let
    post_orders = calculateAllPostOrders rate_mx obs bls site topo
    root_pre_order = PreOrder $
      case topo of
        Leaf i -> getPostOrder $ calculateLeafPostOrder rate_mx obs site i
        Bin {} -> VS.replicate (numOfStates rate_mx) $ 1 / (numOfStates rate_mx)
    total_log_likelihood = log $
      getPreOrder root_pre_order Matrix.<.>
      (getPostOrder $ post_orders ! getNodeId topo)

    go :: MonadState (NodeMap Double) m => PreOrder -> Topology -> m ()
    go pre_order = \case
      Leaf {} -> pure ()
      Bin _ l r -> sequence_ $ do
        let
          subs = tuple2 l r
          sub_pre_orders =
            calculatePreOrder rate_mx bls post_orders pre_order subs
        sub <- subs
        sub_pre_order <- sub_pre_orders
        pure $ do
          let deriv = calculatePartialDerivative
                rate_mx
                total_log_likelihood
                sub_pre_order
                (post_orders ! getNodeId sub)
          modify' $ insert (getNodeId sub) deriv
          go sub_pre_order sub
  in (total_log_likelihood, execState (go root_pre_order topo) mempty)

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
calculatePreOrder rate_mx bls post_orders (PreOrder parent_preorder) subs =
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
      PostOrder sibling_post_order = post_orders ! (getNodeId sibling)
    pure . PreOrder $
      Matrix.tr this_transition_probs Matrix.#>
      (parent_preorder *
        (sibling_transition_probs Matrix.#> sibling_post_order))

-- | Calculate the partial derivative of the log-likelihood w.r.t. the
-- branch length of a given node
calculatePartialDerivative
  :: RateMatrix
  -> Double -- ^ the log-likelihood of the whole tree
  -> PreOrder -- ^ the pre-order traversal of the node
  -> PostOrder -- ^ the post-order traversal of the node
  -> Double -- ^ the partial derivative of the log-likelihood
calculatePartialDerivative (RateMatrix capitalQ) total_ll (PreOrder q) (PostOrder p) =
  (q Matrix.<.> (capitalQ Matrix.#> p)) / exp total_ll

calculateLeafPostOrder
  :: RateMatrix
  -> Observations
  -> Int -- ^ the index of the site
  -> NodeId -- ^ the id of the leaf
  -> PostOrder
calculateLeafPostOrder rate_mx obs site leaf_id =
  let ch = (characters obs ! leaf_id) `VU.unsafeIndex` site
  in PostOrder $ VS.generate (numOfStates rate_mx) $
      \i -> if i == fromIntegral ch then 1 else 0

-- | Calculate all post-order traversals for a given site
calculateAllPostOrders
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Int -- ^ the index of the site
  -> Topology
    -- ^ the tree
  -> NodeMap PostOrder
calculateAllPostOrders rate_mx obs bls site = flip execState mempty . go
  where
    go :: MonadState (NodeMap PostOrder) m => Topology -> m PostOrder
    go topo = do
      po <- case topo of
        Leaf leaf_id -> do
          pure $ calculateLeafPostOrder rate_mx obs site leaf_id
        Bin _ sub1 sub2 -> do
          let
            subs :: Tuple2 Topology
            subs = tuple2 sub1 sub2
          sub_post_orders <- traverse go subs
          pure $ calculatePostOrder rate_mx bls subs sub_post_orders
      modify' $ insert (getNodeId topo) po
      return po
