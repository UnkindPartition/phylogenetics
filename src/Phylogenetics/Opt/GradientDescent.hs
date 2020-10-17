module Phylogenetics.Opt.GradientDescent where

import Phylogenetics.Types
import Phylogenetics.Likelihood_v2

step
  :: Double -- ^ learning rate
  -> Problem
  -> BranchLengths
  -> Maybe (BranchLengths, Double) -- ^ new branch lengths and log-likelihood
step rate prob bls =
  let
    (ll, grad) = gradient prob bls
    bls' = bls + (BranchLength . (rate *) <$> grad)
    ll' = logLikelihood prob bls'
  in
    if
      | ll' > ll -> Just (bls', ll')
      | rate < 1e-8 -> Nothing
      | otherwise -> step (rate / 3) prob bls 
