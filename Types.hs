module Types where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word

newtype NodeId = NodeId Int
  deriving newtype Show

newtype BranchLength = BranchLength Double
  deriving newtype (Eq, Ord, Num, Real, Fractional, RealFrac, Show)

-- | A phylogenetic tree topology.
data Topology
  = Leaf !NodeId
  | Bin !NodeId Topology Topology
  deriving Show

-- | The branch lengths of a phylogenetic tree.
--
-- The map is from the 'NodeId' to the length of a branch leading to that
-- 'NodeId'.
newtype BranchLengths = BranchLengths (IntMap.IntMap BranchLength)
  deriving newtype Show

data Observations = Observations
  { numOfCharacters :: !Int
  , characters :: !(IntMap.IntMap (VU.Vector Word8))
  }
  deriving Show
