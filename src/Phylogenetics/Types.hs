module Phylogenetics.Types where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word
import Numeric.LinearAlgebra

-- | A node in a phylogenetic tree
newtype NodeId = NodeId Int
  deriving newtype Show

-- | The length of a branch in a phylogenetic tree
newtype BranchLength = BranchLength Double
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

-- | The branch lengths of a phylogenetic tree.
--
-- The map is from the 'NodeId' to the length of a branch leading to that
-- 'NodeId'.
newtype BranchLengths = BranchLengths (IntMap.IntMap BranchLength)
  deriving newtype Show

-- | Extract the branch length leading to a given node
getBranchLength :: BranchLengths -> NodeId -> BranchLength
getBranchLength (BranchLengths bl) (NodeId node_id) = bl IntMap.! node_id

-- | A set of observed characters per site. Character states are encoded by
-- integers from 0 to @'numOfStates' - 1@.
data Observations = Observations
  { numOfSites :: !Int
  , characters :: !(IntMap.IntMap (VU.Vector Word8))
  }
  deriving Show

-- | The rate matrix
data RateMatrix = RateMatrix (Matrix Double)

-- | Return the number of states of a character (e.g. 4 for DNA)
numOfStates :: RateMatrix -> Int
numOfStates (RateMatrix mx) = rows mx

-- | Calculate the transition probabilities from one nucleotide to another
transitionProbabilities
  :: RateMatrix
  -> BranchLength
  -> Matrix Double
transitionProbabilities (RateMatrix q) (BranchLength t) = expm (scale t q)
