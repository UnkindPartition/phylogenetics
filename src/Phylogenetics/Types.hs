module Phylogenetics.Types where

import qualified Data.IntMap.Strict as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word
import Control.Monad.State
import Numeric.LinearAlgebra hiding ((<>), fromList)
import Data.Coerce
import GHC.Stack
import GHC.Exts as Exts

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

-- | Return a set of all leaf ids in a topology
leaves :: Topology -> NodeMap ()
leaves = \case
  Leaf leaf_id -> singleton leaf_id ()
  Bin _ l r -> leaves l <> leaves r

allIds :: Topology -> NodeMap ()
allIds = \case
  Leaf i -> singleton i ()
  Bin i l r -> singleton i () <> allIds l <> allIds r

addIdsToTopology :: Topology -> Topology
addIdsToTopology topo0 = evalState (go topo0) 0
  where
    go = \case
      Leaf _ -> Leaf <$> next
      Bin _ l r -> Bin <$> next <*> go l <*> go r
    next = do
      i <- get
      put $! i+1
      return (NodeId i)

-- | The branch lengths of a phylogenetic tree.
--
-- The map is from the 'NodeId' to the length of a branch leading to that
-- 'NodeId'.
type BranchLengths = NodeMap BranchLength

-- | A set of observed characters per site. Character states are encoded by
-- integers from 0 to @'numOfStates' - 1@.
data Observations = Observations
  { numOfSites :: !Int
  , characters :: !(NodeMap (VU.Vector Word8))
  }
  deriving Show

-- | The rate matrix
data RateMatrix = RateMatrix (Matrix Double)
  deriving Show

-- | Return the number of states of a character (e.g. 4 for DNA)
numOfStates :: Num a => RateMatrix -> a
numOfStates (RateMatrix mx) = fromIntegral $ rows mx

-- | Calculate the transition probabilities from one nucleotide to another
transitionProbabilities
  :: RateMatrix
  -> BranchLength
  -> Matrix Double
transitionProbabilities (RateMatrix q) (BranchLength t) = expm (scale t q)

----------------------------------------------------------------------
--                              NodeMap
----------------------------------------------------------------------

newtype NodeMap a = NodeMap (IntMap.IntMap a)
  deriving (Semigroup, Monoid, Show, Functor, Foldable, Traversable)

instance IsList (NodeMap a) where
  type Item (NodeMap a) = (NodeId, a)
  fromList = coerce (fromList @(IntMap.IntMap a))
  toList = coerce (Exts.toList @(IntMap.IntMap a))

lookup :: forall a . NodeId -> NodeMap a -> Maybe a
lookup = coerce (IntMap.lookup @a)

(!) :: forall a . HasCallStack => NodeMap a -> NodeId -> a
(!) = coerce ((IntMap.!) @a)

insert :: forall a . NodeId -> a -> NodeMap a -> NodeMap a
insert = coerce (IntMap.insert @a)

singleton :: forall a . NodeId -> a -> NodeMap a
singleton = coerce (IntMap.singleton @a)

member :: forall a . NodeId -> NodeMap a -> Bool
member = coerce (IntMap.member @a)

size :: forall a . NodeMap a -> Int
size = coerce (IntMap.size @a)

instance Num a => Num (NodeMap a) where
  (+) = coerce (IntMap.unionWith @a (+))
  (-) = coerce (IntMap.unionWith @a (-))
  (*) = coerce (IntMap.unionWith @a (*))
  abs = coerce (IntMap.map @a abs)
  signum = coerce (IntMap.map @a signum)
  negate = coerce (IntMap.map @a negate)
  fromInteger = error "fromInteger is not supported for NodeMaps"
