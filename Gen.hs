-- | Random topology and branch lengths generation
module Gen where

import Data.Random
import Data.Random.Distribution.Binomial
import Data.Random.Distribution.Pareto
import Data.Traversable
import qualified Data.IntMap as IntMap
import Control.Monad.State

import Types

randomTopology
  :: Int -- ^ number of leaves to have in the tree
  -> RVar Topology
randomTopology = flip evalStateT (NodeId 0) . go
  where
    go :: Int -> StateT NodeId RVar Topology
    go n_leaves
      | n_leaves == 1 = Leaf <$> getId
      | n_leaves >= 2 = do
          -- The below generation scheme ensures that both n_leaves_left and
          -- n_leaves_right are positive.
          n_leaves_left <- lift $ (+1) <$> binomial (n_leaves-2) (0.5 :: Double)
          let n_leaves_right = n_leaves - n_leaves_left
          Bin
            <$> getId
            <*> go n_leaves_left
            <*> go n_leaves_right
      | otherwise = undefined

    getId :: MonadState NodeId m => m NodeId
    getId = do
      r@(NodeId i) <- get
      put $! NodeId $! i+1
      return r

randomBranchLength :: RVar BranchLength
randomBranchLength = BranchLength <$> pareto 1 0.9 -- x_min, alpha

randomBranchLengths
  :: Int -- ^ number of leaves in the tree
  -> RVar BranchLengths
randomBranchLengths n_leaves = do
  -- The number of internal nodes in the tree is equal to n_leaves - 1.
  -- The total number of nodes in the tree is 2*n_leaves-1.
  -- The tree root (NodeId 0) does not have an incoming branch, so the
  -- branch ids vary from 1 to 2*n_leaves-2.
  BranchLengths . IntMap.fromList <$>
    for [1 .. 2*n_leaves-2] (\i -> (i,) <$> randomBranchLength)
