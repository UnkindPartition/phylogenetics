module Phylogenetics.Types where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word
import Numeric.LinearAlgebra

-- | A node in a phylogenetic tree
newtype NodeId = NodeId Int
  deriving newtype Show

-- | The time corresponding to a tree branch.
--
-- This is not the «branch length» (the expected number of substitutions),
-- but is linearly related to it for a given rate matrix.
newtype BranchTime = BranchTime Double
  deriving newtype (Eq, Ord, Num, Real, Fractional, RealFrac, Show)

-- | A phylogenetic tree topology.
data Topology
  = Leaf !NodeId
  | Bin !NodeId Topology Topology
  deriving Show

getNodeId :: Topology -> NodeId
getNodeId = \case
  Leaf i -> i
  Bin i _ _ -> i

-- | The time lengths of a phylogenetic tree.
--
-- The map is from the 'NodeId' to the length of a branch leading to that
-- 'NodeId'.
newtype BranchTimes = BranchTimes (IntMap.IntMap BranchTime)
  deriving newtype Show

-- | Extract the branch length leading to a given node
getBranchTime :: BranchTimes -> NodeId -> BranchTime
getBranchTime (BranchTimes bl) (NodeId node_id) = bl IntMap.! node_id

-- | A set of observed characters per site
data Observations = Observations
  { numOfSites :: !Int
  , characters :: !(IntMap.IntMap (VU.Vector Word8))
  }
  deriving Show

-- | The rate matrix
data RateMatrix = RateMatrix (Matrix Double)

-- | Calculate the transition probabilities from one nucleotide to another
transitionProbabilities
  :: RateMatrix
  -> BranchTime
  -> Matrix Double
transitionProbabilities (RateMatrix q) (BranchTime t) = expm (scale t q)
