module Phylogenetics.Types where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word

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

-- | A set of observed characters per site
data Observations = Observations
  { numOfSites :: !Int
  , characters :: !(IntMap.IntMap (VU.Vector Word8))
  }
  deriving Show

-- | The transition matrix stored in the row-major order
data EvolutionModel = EvolutionModel (VU.Vector Double)

-- | Calculate the transition probability from one nucleotide to another
transitionProbability
  :: EvolutionModel
  -> BranchLength
  -> Word8 -- ^ from nucleotide
  -> Word8 -- ^ to nucleotide
  -> Double
transitionProbability (EvolutionModel mx) (BranchLength len) from to =
  exp $ len * VU.unsafeIndex mx (fromIntegral $ from * 4 + to)
