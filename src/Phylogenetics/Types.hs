module Phylogenetics.Types where

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

getBranchLength :: BranchLengths -> NodeId -> BranchLength
getBranchLength (BranchLengths bl) (NodeId node_id) = bl IntMap.! node_id

data Observations = Observations
  { numOfSites :: !Int
  , characters :: !(IntMap.IntMap (VU.Vector Word8))
  }
  deriving Show

-- | The transition matrix stored in the row-major order
data EvolutionModel = EvolutionModel (VU.Vector Double)

transitionProbability
  :: EvolutionModel
  -> BranchLength
  -> Word8 -- from
  -> Word8 -- to
  -> Double
transitionProbability (EvolutionModel mx) (BranchLength len) from to =
  exp $ len * VU.unsafeIndex mx (fromIntegral $ from * 4 + to)
