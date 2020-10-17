module NaiveLikelihood (naiveLikelihood) where

import qualified Data.Vector.Unboxed as VU
import Data.Word
import Numeric.LinearAlgebra as Matrix
import Control.Monad
import Phylogenetics.Types as Phylo

-- | Add all possible values for missing observations
fillObservations
  :: Word8 -- ^ number of states
  -> Topology
  -> Observations
  -> [Observations]
fillObservations n_states topo obs0 =
  case topo of
    Leaf{} -> [obs0]
    Bin node l r -> do
      obs1 <-
        if node `Phylo.member` characters obs0
          then pure obs0
          else do
            these <- replicateM (numOfSites obs0) [0..n_states-1]
            pure obs0
              { characters = insert node (VU.fromList these) (characters obs0) }
      pure obs1
        >>= fillObservations n_states l
        >>= fillObservations n_states r

naiveLikelihood, fullLikelihood
  :: Problem
  -> BranchLengths
  -> Double
naiveLikelihood (Problem rate_mx obs0 topo) bls =
  let
    full_obs = fillObservations (numOfStates rate_mx) topo obs0
    factor =
      -- unless we have a single leaf, we need to multiply by the
      -- probability of observing the tree root
      case topo of
        Leaf{} -> 1
        Bin{} -> (1 / numOfStates rate_mx) ^ (numOfSites obs0)
  in
    factor * sum
      [ fullLikelihood (Problem rate_mx obs topo) bls
      | obs <- full_obs
      ]
fullLikelihood (Problem rate_mx obs@Observations{..} topo) bls =
  case topo of
    Leaf _ -> 1
    Bin parent l@(getNodeId -> l_id) r@(getNodeId -> r_id) ->
      let
        tp_l = transitionProbabilities rate_mx (bls Phylo.! l_id)
        tp_r = transitionProbabilities rate_mx (bls Phylo.! r_id)
      in product
        [ (tp_l
            Matrix.! fromIntegral ((characters Phylo.! parent)  VU.! i)
            Matrix.! fromIntegral ((characters Phylo.! l_id)    VU.! i)) *
          (tp_r
            Matrix.! fromIntegral ((characters Phylo.! parent)  VU.! i)
            Matrix.! fromIntegral ((characters Phylo.! r_id)    VU.! i))
        | i <- [0 .. numOfSites - 1]
        ]
        * fullLikelihood (Problem rate_mx obs l) bls
        * fullLikelihood (Problem rate_mx obs r) bls
