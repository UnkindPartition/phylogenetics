module Phylogenetics.Optimization where

import Phylogenetics.Types
import Phylogenetics.Likelihood_v2

gradientDescentStep
  :: Double -- ^ learning rate
  -> Problem
  -> BranchLengths
  -> Maybe (BranchLengths, Double) -- ^ new branch lengths and log-likelihood
gradientDescentStep rate0 prob bls = go rate0 where
  (ll, grad) = gradient prob bls
  go rate =
    let
      bls' = bls + (BranchLength . (rate *) <$> grad)
      ll' = logLikelihood prob bls'
    in
      if
        | ll' > ll -> Just (bls', ll')
        | rate < 1e-8 -> Nothing
        | otherwise -> go (rate / 3)
