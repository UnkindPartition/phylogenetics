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
               -> BranchLengths
               -> method_state
               -> Maybe (BranchLengths, Double, Double, method_state)
                 -- new bls, actual learning rate, new likelihood, new state
  , methodInit :: method_state
  }

gradientDescent
  :: Double -- ^ learning rate
  -> Method
gradientDescent rate = Method
  { methodName = printf "Gradient descent; rate = %.1g" rate
  , methodStep = gradientDescentStep
  , methodInit = rate
  }

gradientDescentStep
  :: Problem
  -> BranchLengths
  -> Double -- ^ learning rate
  -> Maybe (BranchLengths, Double, Double, Double) -- ^ new branch lengths, actual learning rate, log-likelihood, new state
gradientDescentStep prob bls rate0 = go rate0 where
  (ll, grad) = gradient prob bls
  go rate =
    let
      bls' = (max 1e-10) <$> bls + (BranchLength . (rate *) <$> grad)
      ll' = logLikelihood prob bls'
    in
      if
        | ll' > ll -> Just (bls', rate, ll', rate)
        | rate < 1e-10 -> Nothing
        | otherwise -> go (rate / 3)
