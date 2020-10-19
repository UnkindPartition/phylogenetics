module Phylogenetics.Optimization
  ( Method(..)
  , gradientDescent
  , noDescent
  , noDescentAccelerated
  )
  where

import Text.Printf (printf)
import Phylogenetics.Types
import Phylogenetics.Likelihood_v2

data Method = forall method_state . Method
  { methodName :: String
  , methodStep :: Problem
               -> method_state
               -> Maybe (BranchLengths, Double, Double, Double, method_state)
                 -- new bls, actual learning rate, new likelihood, gradient norm, new state
  , methodInit :: Problem -> BranchLengths -> method_state
    -- ^ construct the initial state from the starting point
  }

gradientUpdate
  :: Double -- ^ learning rate
  -> BranchLengths
  -> Gradient
  -> BranchLengths
gradientUpdate rate bls grad =
  max 1e-10 <$> (bls + (BranchLength . (rate *) <$> grad))

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
  -> Maybe (BranchLengths, Double, Double, Double, BranchLengths)
gradientDescentStep rate0 prob bls = go rate0 where
  go rate =
    let
      (ll, grad) = gradient prob bls
      bls' = gradientUpdate rate bls grad
      ll' = logLikelihood prob bls'
    in
      if
        | ll' > ll -> Just (bls', rate, ll', l2norm grad, bls')
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
  { methodName = "Without descent"
  , methodInit = \prob bls -> NoDescentState 1e-4 (1/0) (bls, snd $ gradient prob bls) Nothing
  , methodStep = \prob (NoDescentState prev_λ prev_θ (this_bls, this_grad) mb_prev) ->
      case mb_prev of
        Nothing -> -- first step
          let
            new_bls = gradientUpdate prev_λ this_bls this_grad
            (new_ll, new_grad) = gradient prob new_bls
          in
            Just ( new_bls, prev_λ, new_ll, l2norm this_grad
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
            Just ( new_bls, λ, new_ll, l2norm this_grad
                 , NoDescentState λ θ (new_bls, new_grad) (Just (this_bls, this_grad)))
  }

data NoDescentAccState = NoDescentAccState
  Double -- little λ
  Double -- little θ
  Double -- big Λ
  Double -- big Θ
  BranchLengths -- y
  (BranchLengths, Gradient)
  (Maybe (BranchLengths, Gradient))
    -- previous point and gradient

noDescentAccelerated :: Method
noDescentAccelerated = Method
  { methodName = "Adaptive accelerated no-descent"
  , methodInit = \prob bls -> NoDescentAccState 1e-4 1e-4 (1/0) (1/0) bls (bls, snd $ gradient prob bls) Nothing
  , methodStep = \prob (NoDescentAccState this_λ this_θ this_Λ this_Θ this_y (this_x, this_grad) mb_prev) ->
      case mb_prev of
        Nothing -> -- first step
          let
            new_x = gradientUpdate this_λ this_x this_grad
            (new_ll, new_grad) = gradient prob new_x
          in
            Just ( new_x, this_λ, new_ll, l2norm this_grad
                 , NoDescentAccState this_λ this_θ this_Λ this_Θ new_x (new_x, new_grad) (Just (this_x, this_grad)))
        Just (prev_x, prev_grad) ->
          let
            new_λ = min
              (sqrt (1+this_θ) * this_λ)
              (l2norm (this_x - prev_x) / (2 * l2norm (this_grad - prev_grad)))
            new_Λ = min
              (sqrt (1+this_Θ) * this_Λ)
              (l2norm (this_grad - prev_grad) / (2 * l2norm (this_x - prev_x)))

            inv_sqrt_λ = 1 / sqrt new_λ
            sqrt_Λ = sqrt new_Λ
            β = (inv_sqrt_λ - sqrt_Λ) / (inv_sqrt_λ + sqrt_Λ)

            new_y = gradientUpdate new_λ this_x this_grad
            new_x = max 1e-10 <$> (new_y + ((BranchLength β *) <$> (new_y - this_y)))

            new_θ = new_λ / this_λ
            new_Θ = new_Λ / this_Λ

            (new_ll, new_grad) = gradient prob new_x
          in
            Just ( new_x, new_λ, new_ll, l2norm this_grad
                 , NoDescentAccState new_λ new_θ new_Λ new_Θ new_y (new_x, new_grad) (Just (this_x, this_grad)))
  }
