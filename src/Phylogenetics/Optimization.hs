module Phylogenetics.Optimization
  ( Method(..)
  , gradientDescent
  )
  where

import Text.Printf (printf)
import Phylogenetics.Types
import Phylogenetics.Likelihood_v2

data Method = forall method_state . Method
  { methodName :: String
  , methodStep :: Problem
               -> method_state
               -> Maybe (BranchLengths, Double, Double, method_state)
                 -- new bls, actual learning rate, new likelihood, new state
  , methodInit :: BranchLengths -> method_state
    -- ^ construct the initial state from the starting point
  }

gradientUpdate
  :: Double -- ^ learning rate
  -> BranchLengths
  -> Gradient
  -> BranchLengths
gradientUpdate rate bls grad =
  (max 1e-10) <$> bls + (BranchLength . (rate *) <$> grad)

gradientDescent
  :: Double -- ^ learning rate
  -> Method
gradientDescent rate = Method
  { methodName = printf "Gradient descent; rate = %.1g" rate
  , methodStep = gradientDescentStep rate
  , methodInit = id
  }

gradientDescentStep
  :: Double -- ^ learning rate
  -> Problem
  -> BranchLengths
  -> Maybe (BranchLengths, Double, Double, BranchLengths) -- ^ new branch lengths, actual learning rate, log-likelihood, new state
gradientDescentStep rate0 prob bls = go rate0 where
  go rate =
    let
      (ll, grad) = gradient prob bls
      bls' = gradientUpdate rate bls grad
      ll' = logLikelihood prob bls'
    in
      if
        | ll' > ll -> Just (bls', rate, ll', bls')
        | rate < 1e-10 -> Nothing
        | otherwise -> go (rate / 5)
