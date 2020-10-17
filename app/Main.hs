import Options.Applicative
import Control.Monad
import Control.Monad.State
import Control.Monad.Trans.Maybe
import Control.DeepSeq
import Control.Exception
import Data.Random
import Text.Printf (printf)
import System.IO
import System.Clock
import Phylogenetics.Types
import Phylogenetics.Gen
import Phylogenetics.Likelihood_v2
import Phylogenetics.Optimization

main :: IO ()
main = join . customExecParser (prefs showHelpOnError) $
  info (helper <*> parser)
  (  fullDesc
  )
  where
    parser :: Parser (IO ())
    parser =
      hsubparser
        ( command "trace" (info traceParser fullDesc) )

traceParser :: Parser (IO ())
traceParser = trace
  <$> option auto
    (  long "num-steps"
    <> metavar "NUMBER"
    <> help "number of optimization steps"
    <> value 30
    <> showDefault
    )

data TraceState method_state = TraceState
  { trace_method_state :: !method_state
  }

trace :: Int -> IO ()
trace num_steps = do
  hSetBuffering stdout LineBuffering

  (prob, true_bls) <- sample $ do
    rate_mx <- rateMatrix bd
    topo <- topology bd
    bls <- branchLengths bd topo
    obs <- realisticObservations bd rate_mx topo bls
    return (Problem rate_mx obs topo, bls)

  let init_bls = 0.1 <$ true_bls

  printf "method,step,rate,ll,mse,time\n"

  forM_ methods $ \Method{..} -> do
    printf "%s,0,NA,%.6f,%.6g,0\n"
      methodName
      (logLikelihood prob init_bls)
      (calculateMSE init_bls true_bls)
    flip evalStateT (TraceState (methodInit prob init_bls)) $
      forM_ [1 .. num_steps] $ \istep ->
      runMaybeT $ do
        prev_state <- get
        t0 <- liftIO $ getTime ProcessCPUTime
        let step_result = methodStep prob (trace_method_state prev_state)
        case step_result of
          Nothing -> pure ()
          Just (bls, actual_learning_rate, ll, _) -> do
            liftIO . evaluate . rnf $ bls
            liftIO . evaluate . rnf $ actual_learning_rate
            liftIO . evaluate . rnf $ ll
        t1 <- liftIO $ getTime ProcessCPUTime
        let
          time :: Double
          time = fromIntegral (toNanoSecs (diffTimeSpec t1 t0)) / 1e9
        case methodStep prob (trace_method_state prev_state) of
          Just (bls, actual_learning_rate, ll, st) -> do
            put $! TraceState st
            liftIO $ printf "%s,%d,%.4g,%.6f,%.6g,%.3g\n" methodName istep actual_learning_rate ll (calculateMSE bls true_bls) time
            --liftIO $ print bls
          Nothing -> mzero

  return ()
  where
    bd = dnaBaseDistributions
      { numberOfSitesDistribution = pure 300
      , numberOfLeavesInTreeDistribution = pure 30
      }
    methods =
      [noDescent] ++ (gradientDescent <$> [1e-3, 1e-4])

calculateMSE
  :: BranchLengths
  -> BranchLengths
  -> Double
calculateMSE bls true_bls =
  case sum . fmap (^(2::Int)) $ bls - true_bls of
    BranchLength sum_errors_squared ->
      sum_errors_squared / fromIntegral (size true_bls)
