{-# LANGUAGE OverloadedLists #-}
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck as QC hiding ((><), scale)
import Test.Tasty.ExpectedFailure
import Text.Printf (printf)
import Data.Tuple.Homogenous
import Numeric.LinearAlgebra as Matrix
import Phylogenetics.Types as Phylo
import qualified Phylogenetics.Likelihood_v1 as V1
import qualified Phylogenetics.Likelihood_v2 as V2
import qualified Phylogenetics.Gen as Gen
import qualified Phylogenetics.Optimization as Opt
import NaiveLikelihood
import NaiveGradient

qcDistributions :: Gen.BaseDistributions Gen
qcDistributions = Gen.BaseDistributions
  { numberOfLeavesInTreeDistribution = frequency $
      (5,pure 1) : [(i^(2::Integer), pure $ 6-i) | i <- [1..4]]
  , numberOfLeavesInLeftSubtreeDistribution = \n -> choose (1, n-1)
  , numberOfSitesDistribution = choose (1, 3)
  , numberOfCharacterStatesDistribution = choose (2, 4)
  , characterUniformDistribution = \n -> choose (0, n-1)
  , characterCategoricalDistribution = \freqs -> frequency $
      zip (round . (1e6*) <$> freqs) (pure <$> [0..])
  , branchLengthDistribution = choose (0,0.5)
  , rateDistribution = getPositive <$> arbitrary
  }

instance Arbitrary RateMatrix where
  arbitrary = Gen.rateMatrix qcDistributions
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
  arbitrary = Gen.branchLength qcDistributions

instance Arbitrary Topology where
  arbitrary = Gen.topology qcDistributions
  shrink = \case
    Leaf{} -> []
    Bin _ l r -> [l, r]

main :: IO ()
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
          bls =
            [ (NodeId 0, 3.0)
            , (NodeId 2, 5.0)
            ]
          obs = Observations 1
            [ (NodeId 0, [1])
            , (NodeId 2, [0])
            ]
          prob = Problem (RateMatrix rm) obs topo
          ll = V1.logLikelihood prob bls
          tp0 = expm (scale 3.0 rm)
          tp2 = expm (scale 5.0 rm)
          expected_ll = log $
            0.5 * (tp0 Matrix.! 0 Matrix.! 1) * (tp2 Matrix.! 0 Matrix.! 0) +
            0.5 * (tp0 Matrix.! 1 Matrix.! 1) * (tp2 Matrix.! 1 Matrix.! 0)
        assertBool
          (printf "Expected %.3f, got %.3f" expected_ll ll)
          (abs (expected_ll - ll) < 1e-10)
    , testProperty "naiveLikelihood vs V1.logLikelihood" $
        testLikelihoodCalculation $ Tuple2
          ( \prob bls -> log $ naiveLikelihood prob bls
          , V1.logLikelihood
          )
    , testProperty "V1.logLikelihood vs V2.logLikelihood" $
        testLikelihoodCalculation $ Tuple2
          ( V1.logLikelihood
          , V2.logLikelihood
          )
    , testProperty "V2.logLikelihood vs V2.gradient" $
        testLikelihoodCalculation $ Tuple2
          ( V2.logLikelihood
          , \prob bls -> fst $ V2.gradient prob bls
          )
    , testProperty "naiveGradient vs V2.gradient" $
        testWithPhylogeneticsModel $ \prob bls ->
          let
            naive_gr = naiveGradient prob bls
            gr = snd $ V2.gradient prob bls
            l1_diff = sum (abs (naive_gr - gr))
          in
            counterexample (show l1_diff) $
            counterexample (show naive_gr) $
            counterexample (show gr) $
              l1_diff < 1e-4
    ]
  , testGroup "Optimization"
    [ testGroup "GradientDescent"
      [ testProperty "One step improves the likelihood of the original branch lengths" $
          testWithRealisticPhylogeneticsModel $ \prob bls ->
            case Opt.gradientDescentStep 0.01 prob bls of
              Nothing -> label "step returned Nothing" True
              Just (bls', ll) -> property $
                (ll - V2.logLikelihood prob bls') < 1e-10 &&
                (ll > V2.logLikelihood prob bls)
      , testProperty "One step improves the likelihood of initial branch lengths" $
          testWithRealisticPhylogeneticsModel $ \prob bls ->
            let
              init_bls = 0.1 <$ bls
              init_ll = V2.logLikelihood prob init_bls
              orig_ll = V2.logLikelihood prob bls
            in
            counterexample (printf "init_ll = %.3f, orig_ll = %.3f" init_ll orig_ll) $
            case Opt.gradientDescentStep 0.01 prob init_bls of
              Nothing -> counterexample "step returned Nothing" False
              Just (_bls', step_ll) -> property $
                  counterexample (printf "step_ll = %.3f" step_ll) $
                  classify (step_ll > orig_ll) "Improved above the original LL" $
                    init_ll < step_ll
      ]
    ]
  ]

testWithRealisticPhylogeneticsModel
  :: QC.Testable a
  => (Problem -> BranchLengths -> a)
  -> Property
testWithRealisticPhylogeneticsModel k =
  forAll (Gen.rateMatrix bd) $ \rate_mx ->
  forAll (Gen.topology bd) $ \topo ->
  forAll (Gen.branchLengths bd topo) $ \bls ->
  forAll (Gen.realisticObservations bd rate_mx topo bls) $ \obs ->
  k (Problem rate_mx obs topo) bls
  where
    bd = qcDistributions
      { Gen.numberOfSitesDistribution = pure 100
      , Gen.numberOfLeavesInTreeDistribution = choose (3,10)
      }

testWithPhylogeneticsModel
  :: QC.Testable a
  => (Problem -> BranchLengths -> a)
  -> Property
testWithPhylogeneticsModel k =
  forAll (Gen.problem qcDistributions) $ \prob@(Problem _ _ topo) ->
  forAll (Gen.branchLengths qcDistributions topo) $ \bls ->
  k prob bls

testLikelihoodCalculation
  :: Tuple2 (Problem -> BranchLengths -> Double)
     -- ^ the two log-likelihood calculation functins
  -> Property
testLikelihoodCalculation ll_fns = testWithPhylogeneticsModel $ \prob@(Problem _ _ topo) bls ->
  let
    Tuple2 (ll1, ll2) = do
      ll_fn <- ll_fns
      pure $ ll_fn prob bls
  in
    label (show (Phylo.size $ leaves topo) ++ " leaves") $
    counterexample (printf "First LL: %.3f, second LL: %.3f" ll1 ll2) $
      (abs (ll1 - ll2) < 1e-10)

isStochastic :: Matrix Double -> Bool
isStochastic mx =
  let row_sums = mx #> konst 1 (cols mx)
  in norm_Inf (row_sums - konst 1 (rows mx)) < 1e-10
