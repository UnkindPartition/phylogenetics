module Phylogenetics.Optimization
  ( Method(..)
  , gradientDescent
  , noDescent
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
  , methodInit :: Problem -> BranchLengths -> method_state
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
  , methodInit = \_prob bls -> bls
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

data NoDescentState = NoDescentState
  Double -- λ
  Double -- θ
  (BranchLengths, Gradient)
  (Maybe (BranchLengths, Gradient))
    -- previous point and gradient

noDescent :: Method
noDescent = Method
  { methodName = printf "Without descent"
  , methodInit = \prob bls -> NoDescentState 1e-10 (1/0) (bls, snd $ gradient prob bls) Nothing
  , methodStep = \prob (NoDescentState prev_λ prev_θ (this_bls, this_grad) mb_prev) ->
      case mb_prev of
        Nothing -> -- first step
          let
            new_bls = gradientUpdate prev_λ this_bls this_grad
            (new_ll, new_grad) = gradient prob new_bls
          in
            Just ( new_bls, prev_λ, new_ll
                 , NoDescentState prev_λ prev_θ (new_bls, new_grad) (Just (this_bls, this_grad)))
        Just (prev_bls, prev_grad) ->
          let
            new_bls = gradientUpdate λ this_bls this_grad
            λ = min
              (sqrt (1+prev_θ) * prev_λ)
              (l2norm (this_bls - prev_bls) / (2 * l2norm (this_grad - prev_grad)))
            (new_ll, new_grad) = gradient prob new_bls
            θ = λ / prev_λ
          in
            Just ( new_bls, λ, new_ll
                 , NoDescentState λ θ (new_bls, new_grad) (Just (this_bls, this_grad)))
  }
