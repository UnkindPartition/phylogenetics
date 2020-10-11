{-# LANGUAGE OverloadedLists #-}
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck hiding ((><), scale)
import Test.Tasty.ExpectedFailure
import Test.QuickCheck.Gen hiding (scale)
import Text.Printf (printf)
import qualified Data.IntSet as IntSet
import qualified Data.IntMap as IntMap
import qualified Data.Vector.Unboxed as VU
import Data.Word
import Control.Monad
import Numeric.LinearAlgebra
import Phylogenetics.Types
import Phylogenetics.Likelihood
import NaiveLikelihood

rateMatrixGen
  :: Int -- ^ size, i.e. the number of states
  -> Gen RateMatrix
rateMatrixGen size = do
  rows <- forM [0 .. size - 1] $ \i -> do
    rates0 <- replicateM (size - 1) (getPositive <$> arbitrary)
    let rates = map (/ sum rates0) rates0
    return $ fromList $ take i rates ++ [- 1] ++ drop i rates
  return $ RateMatrix $ fromRows rows

instance Arbitrary RateMatrix where
  arbitrary = do
    size <- choose (2, 30)
    rateMatrixGen size
  shrink (RateMatrix mx)
    | old_size > 2 =
        let
          size = if old_size > 5 then old_size `div` 2 else old_size - 1
          rows = (\f -> zipWith f [0 .. size - 1] (toRows mx)) $ \i old_row ->
            let
              rates0 = take (size - 1) $ filter (> 0) $ toList old_row
              rates = map (/ sum rates0) rates0
            in
              fromList $ take i rates ++ [- 1] ++ drop i rates
        in return $ RateMatrix $ fromRows rows
    | otherwise = []
    where old_size = cols mx

instance Arbitrary BranchLength where
  arbitrary = BranchLength . getPositive <$> arbitrary

-- | Generate branch lengths for a given topology
branchLengthsGen
  :: Topology
  -> Gen BranchLengths
branchLengthsGen topo = do
  let node_ids = allIds topo
  fmap BranchLengths . sequence
    . IntMap.fromList . map (, arbitrary)
    . IntSet.toList $ node_ids

instance Arbitrary Topology where
  arbitrary = addIdsToTopology <$> frequency
    [ (1, Bin (NodeId 0) <$> arbitrary <*> arbitrary)
    , (5, pure $ Leaf $ NodeId 0)
    ]

observationsGen :: RateMatrix -> Topology -> Gen Observations
observationsGen rate_mx topo = do
  let
    leaf_ids = leaves topo
    m = numOfStates rate_mx

  n <- choose (1, 3) -- number of sites

  let
    characterAtSite :: Gen (VU.Vector Word8)
    characterAtSite = VU.fromList <$> vectorOf n (choose (0, m-1))

  fmap (Observations n) . sequence
    . IntMap.fromList . map (, characterAtSite)
    . IntSet.toList $ leaf_ids


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
        assertBool
          (printf "Expected %.3f, got %.3f" expected_ll ll)
          (abs (expected_ll - ll) < 1e-10)
    , testProperty "Checking against naive implmenetation" $ \topo rate_mx ->
        forAll (branchLengthsGen topo) $ \bls ->
        forAll (observationsGen rate_mx topo) $ \obs ->
          let
            naive_ll = log $ naiveLikelihood rate_mx obs bls topo
            ll       = logLikelihood         rate_mx obs bls topo
          in
            counterexample (printf "Naive LL: %.3f, efficient: %.3f" naive_ll ll) $
              (abs (naive_ll - ll) < 1e-10)
    ]
  ]

isStochastic :: Matrix Double -> Bool
isStochastic mx =
  let row_sums = mx #> konst 1 (cols mx)
  in norm_Inf (row_sums - konst 1 (rows mx)) < 1e-10
