module NaiveLikelihood (naiveLikelihood) where

import qualified Data.IntMap as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word
import Numeric.LinearAlgebra
import Control.Monad
import qualified Control.Foldl as L
import Phylogenetics.Types

-- | Add all possible values for missing observations
fillObservations
  :: Word8 -- ^ number of states
  -> Topology
  -> Observations
  -> [Observations]
fillObservations n_states topo obs0 =
  case topo of
    Leaf{} -> [obs0]
    Bin (NodeId node) l r -> do
      obs1 <-
        if node `IntMap.member` characters obs0
          then pure obs0
          else do
            these <- replicateM (numOfSites obs0) [0..n_states-1]
            pure obs0
              { characters = IntMap.insert node (VU.fromList these) (characters obs0) }
      pure obs1
        >>= fillObservations n_states l
        >>= fillObservations n_states r

naiveLikelihood, fullLikelihood
  :: RateMatrix
  -> Observations
  -> BranchLengths
  -> Topology
  -> Double
naiveLikelihood rate_mx obs0 bls topo =
  let
    full_obs = fillObservations (numOfStates rate_mx) topo obs0
    average = (/) <$> L.sum <*> L.genericLength
  in
    L.fold average
      [ fullLikelihood rate_mx obs bls topo
      | obs <- full_obs
      ]
fullLikelihood rate_mx obs@Observations{..} bls =
  \case
    Leaf _ -> 1
    Bin (NodeId parent) l@(getNodeId -> NodeId l_id) r@(getNodeId -> NodeId r_id) ->
      let
        tp_l = transitionProbabilities rate_mx (getBranchLength bls (NodeId l_id))
        tp_r = transitionProbabilities rate_mx (getBranchLength bls (NodeId r_id))
      in product
        [ (tp_l
            ! fromIntegral ((characters IntMap.! parent)  VU.! i)
            ! fromIntegral ((characters IntMap.! l_id)    VU.! i)) *
          (tp_r
            ! fromIntegral ((characters IntMap.! parent)  VU.! i)
            ! fromIntegral ((characters IntMap.! r_id)    VU.! i))
        | i <- [0 .. numOfSites - 1]
        ]
        * fullLikelihood rate_mx obs bls l
        * fullLikelihood rate_mx obs bls r
