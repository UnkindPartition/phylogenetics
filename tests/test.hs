{-# LANGUAGE OverloadedLists #-}
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck hiding ((><), scale)
import Test.Tasty.ExpectedFailure
import Test.QuickCheck.Gen hiding (scale)
import Text.Printf (printf)
import Control.Monad
import Numeric.LinearAlgebra
import Phylogenetics.Types
import Phylogenetics.Likelihood

instance Arbitrary RateMatrix where
  arbitrary = do
    size <- choose (2, 30)
    rows <- forM [0 .. size - 1] $ \i -> do
      rates <- replicateM (size - 1) (getPositive <$> arbitrary)
      return $ fromList $ take i rates ++ [- sum rates] ++ drop i rates
    return $ RateMatrix $ fromRows rows

instance Arbitrary BranchLength where
  arbitrary = BranchLength . getPositive <$> arbitrary

main = defaultMain $ testGroup "Tests"
  [ testGroup "transitionProbabilities"
    [ testProperty "transitionProbabilities is a stochastic matrix" $ \rm bl ->
        isStochastic (transitionProbabilities rm bl)
    , ignoreTestBecause "https://github.com/haskell-numerics/hmatrix/issues/341" $
        testProperty "transitionProbabilities at 0 is identity" $ \rm ->
          norm_Inf (transitionProbabilities rm (BranchLength 0) - ident (numOfStates rm)) < 1e-10
    ]
  , testGroup "Likelihood"
    [ testCase "Simple 2-state model" $ do
        let
          rm, tp0, tp2 :: Matrix Double
          rm = (2><2) [-1, 1, 1, -1]
          -- Tree: 1 -> 0, 1 -> 2
          topo = Bin (NodeId 1) (Leaf (NodeId 0)) (Leaf (NodeId 2))
          bls = BranchLengths
            [ (0, 3.0)
            , (2, 5.0)
            ]
          obs = Observations 1
            [ (0, [1])
            , (2, [0])
            ]
          ll = logLikelihood (RateMatrix rm) obs bls topo
          tp0 = expm (scale 3.0 rm)
          tp2 = expm (scale 5.0 rm)
          expected_ll = log $
            0.5 * (tp0 ! 0 ! 1) * (tp2 ! 0 ! 0) +
            0.5 * (tp0 ! 1 ! 1) * (tp2 ! 1 ! 0)
        unless (abs (expected_ll - ll) < 1e-10) $
          assertFailure $ printf "Expected %.3f, got %.3f" expected_ll ll
    ]
  ]

isStochastic :: Matrix Double -> Bool
isStochastic mx =
  let row_sums = mx #> konst 1 (cols mx)
  in norm_Inf (row_sums - konst 1 (rows mx)) < 1e-10
