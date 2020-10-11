-- | Generate random phylogenetic trees and observations.
--
-- The functions here abstract over the monadic RNG to accommodate both
-- random-fu (for numeric experiments) and QuickCheck (for testing).
--
-- The module is intended to be imported qualified.
module Phylogenetics.Gen where

import Data.Word
import qualified Data.IntSet as IntSet
import qualified Data.IntMap as IntMap
import qualified Data.Vector.Unboxed as VU
import Numeric.LinearAlgebra (fromList, fromRows)
import Control.Monad

import Phylogenetics.Types

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
  , characterStateDistribution :: Word8 -> m Word8
    -- ^ Generate a single random character, given the total number of
    -- characters (should be between 0 and n-1)
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
branchLengths bd topo = do
  let node_ids = allIds topo
  fmap BranchLengths . sequence
    . IntMap.fromList . map (, branchLength bd)
    . IntSet.toList $ node_ids

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
    characterAtSite = VU.fromList <$> replicateM n (characterStateDistribution m)

  fmap (Observations n) . sequence
    . IntMap.fromList . map (, characterAtSite)
    . IntSet.toList $ leaf_ids

rateMatrix
  :: forall m . Monad m
  => BaseDistributions m
  -> m RateMatrix
rateMatrix BaseDistributions{..} = do
  size <- fromIntegral <$> numberOfCharacterStatesDistribution
  rows <- forM [0 .. size - 1] $ \i -> do
    rates0 <- replicateM (size - 1) rateDistribution
    let rates = map (/ sum rates0) rates0
    return $ fromList $ take i rates ++ [- 1] ++ drop i rates
  return $ RateMatrix $ fromRows rows
