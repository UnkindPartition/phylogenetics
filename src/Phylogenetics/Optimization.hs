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
  , methodStep = gradientDescentStep rate
  , methodInit = ()
  }

gradientDescentStep
  :: Double -- ^ learning rate
  -> Problem
  -> BranchLengths
  -> () -- ^ no state
  -> Maybe (BranchLengths, Double, Double, ()) -- ^ new branch lengths, actual learning rate, log-likelihood, new state
gradientDescentStep rate0 prob bls () = go rate0 where
  (ll, grad) = gradient prob bls
  go rate =
    let
      bls' = (max 1e-10) <$> bls + (BranchLength . (rate *) <$> grad)
      ll' = logLikelihood prob bls'
    in
      if
        | ll' > ll -> Just (bls', rate, ll', ())
        | rate < 1e-8 -> Nothing
        | otherwise -> go (rate / 3)
