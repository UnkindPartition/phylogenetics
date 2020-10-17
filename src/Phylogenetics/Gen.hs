-- | Generate random phylogenetic trees and observations.
--
-- The functions here abstract over the monadic RNG to accommodate both
-- random-fu (for numeric experiments) and QuickCheck (for testing).
--
-- The module is intended to be imported qualified.
module Phylogenetics.Gen where

import Data.Word
import qualified Data.Vector.Unboxed as VU
import Numeric.LinearAlgebra as Matrix (fromList, toList, fromRows, (!))
import Control.Monad
import Control.Monad.State

import Phylogenetics.Types as Phylo hiding (size)

data BaseDistributions m = BaseDistributions
  { numberOfLeavesInTreeDistribution :: m Int
    -- ^ How many leaves are in a tree?
  , numberOfLeavesInLeftSubtreeDistribution :: Int -> m Int
    -- ^ Given the total number of leaves in a tree,
    -- how many of them are in the left subtree?
    -- The result should be between 1 and n-1.
  , numberOfSitesDistribution :: m Int
    -- ^ How many sites are we observing?
  , numberOfCharacterStatesDistribution :: m Word8
    -- ^ How many states does each character have?
  , characterUniformDistribution :: Word8 -> m Word8
    -- ^ Generate a single random character, given the total number of
    -- characters (should be between 0 and n-1)
  , characterCategoricalDistribution :: [Double] -> m Word8
    -- ^ Generate a single random character from the given categorical
    -- distribution
  , branchLengthDistribution :: m Double
    -- ^ Branch length; should be > 0.
  , rateDistribution :: m Double
    -- ^ Distribution of rates in the rate matrix. Note that the rates will
    -- be normalized such that their sum (excluding the diagonal entry) in
    -- each row is 1. Rates must be positive.
  }

topology
  :: forall m . Monad m
  => BaseDistributions m
  -> m Topology
topology BaseDistributions{..} = do
  n_leaves <- numberOfLeavesInTreeDistribution
  fmap addIdsToTopology $ go n_leaves
  where
    go :: Int -> m Topology
    go n_leaves
      | n_leaves == 1 = pure $ Leaf (NodeId 0)
      | n_leaves >= 2 = do
          n_leaves_left <- numberOfLeavesInLeftSubtreeDistribution n_leaves
          let n_leaves_right = n_leaves - n_leaves_left
          Bin (NodeId 0)
            <$> go n_leaves_left
            <*> go n_leaves_right
      | otherwise = undefined

branchLength
  :: forall m . Monad m
  => BaseDistributions m
  -> m BranchLength
branchLength BaseDistributions{..} = BranchLength <$> branchLengthDistribution

branchLengths
  :: forall m . Monad m
  => BaseDistributions m
  -> Topology
  -> m BranchLengths
branchLengths bd = traverse (const $ branchLength bd) . allIds

-- | Generate a random set of observations that doesn't really take into
-- account the rate matrix or branch lengths
observations
  :: forall m . Monad m
  => BaseDistributions m
  -> RateMatrix
  -> Topology
  -> m Observations
observations BaseDistributions{..} rate_mx topo = do
  let
    leaf_ids = leaves topo
    m = numOfStates rate_mx

  n <- numberOfSitesDistribution

  let
    characterAtSite :: m (VU.Vector Word8)
    characterAtSite = VU.fromList <$> replicateM n (characterUniformDistribution m)

  fmap (Observations n) . traverse (const characterAtSite) $ leaf_ids

-- | Generate a realistic set of observations that takes into account the
-- rate matrix and branch lengths
realisticObservations
  :: forall m . Monad m
  => BaseDistributions m
  -> RateMatrix
  -> Topology
  -> BranchLengths
  -> m Observations
realisticObservations BaseDistributions{..} rate_mx topo0 bls = do
  let
    m = numOfStates rate_mx
  numOfSites <- numberOfSitesDistribution
  charactersAtRoot <- VU.fromList <$> replicateM numOfSites (characterUniformDistribution m)
  characters <- execStateT (go charactersAtRoot topo0) mempty
  pure Observations{..}
  where
    go :: VU.Vector Word8 -> Topology -> StateT (NodeMap (VU.Vector Word8)) m ()
    go chars = \case
        Leaf leaf_id -> modify' $ insert leaf_id chars
        Bin _ l r -> descend chars l >> descend chars r
    descend :: VU.Vector Word8 -> Topology -> StateT (NodeMap (VU.Vector Word8)) m ()
    descend parent_chars topo = do
      let
        bl = bls Phylo.! getNodeId topo
        probs = transitionProbabilities rate_mx bl
      these_chars <- VU.forM parent_chars $ \char ->
        (lift . characterCategoricalDistribution . Matrix.toList)
        (probs Matrix.! fromIntegral char)
      go these_chars topo

rateMatrix
  :: forall m . Monad m
  => BaseDistributions m
  -> m RateMatrix
rateMatrix BaseDistributions{..} = do
  size <- fromIntegral <$> numberOfCharacterStatesDistribution
  rows <- forM [0 .. size - 1] $ \i -> do
    rates0 <- replicateM (size - 1) rateDistribution
    let rates = map (/ sum rates0) rates0
    return $ Matrix.fromList $ take i rates ++ [- 1] ++ drop i rates
  return $ RateMatrix $ fromRows rows

problem
  :: forall m . Monad m
  => BaseDistributions m
  -> m Problem
problem bd = do
  rate_mx <- rateMatrix bd
  topo <- topology bd
  obs <- observations bd rate_mx topo
  pure $ Problem rate_mx obs topo
