import Criterion
import Criterion.Main
import Data.Random
import Data.Random.Distribution.Binomial
import Data.Random.Distribution.Uniform
import Data.Random.Distribution.Exponential
import qualified Phylogenetics.Likelihood_v1 as V1
import qualified Phylogenetics.Likelihood_v2 as V2
import qualified Phylogenetics.Gen as Gen
import Phylogenetics.Types

benchDistributions :: Gen.BaseDistributions RVar
benchDistributions = Gen.BaseDistributions
  { numberOfLeavesInTreeDistribution = pure 100
  , numberOfLeavesInLeftSubtreeDistribution = \n -> (1 +) <$> binomial (n-2) (0.5 :: Double)
  , numberOfSitesDistribution = pure 1
  , numberOfCharacterStatesDistribution = pure 4
  , characterStateDistribution = \n -> integralUniform 0 (n-1)
  , branchLengthDistribution = exponential 0.2
  , rateDistribution = stdUniform
  }

main :: IO ()
main = do
  prob@(Problem _ _ topo) <- sample $ Gen.problem benchDistributions
  bls <- sample $ Gen.branchLengths benchDistributions topo
  defaultMain
    [ bench "V1.logLikelihood" $ whnf (V1.logLikelihood prob) bls
    , bench "V2.logLikelihood" $ whnf (V2.logLikelihood prob) bls
    , bench "V2.gradient" $ whnf (fst . V2.gradient prob) bls
    ]
