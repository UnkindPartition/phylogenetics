import Options.Applicative
import Control.Monad
import Control.Monad.State
import Control.Monad.Trans.Maybe
import Control.DeepSeq
import Control.Exception
import Data.Random
import Data.Random.Source.PureMT
import Data.IORef
import Data.Word
import Text.Printf
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
  <*> option auto
    (  long "seed"
    <> metavar "NUMBER"
    <> help "random seed"
    <> value 2020
    <> showDefault
    )
  <*> option auto
    (  long "num-sites"
    <> metavar "NUMBER"
    <> help "Number of observed sites"
    <> value 100
    <> showDefault
    )
  <*> option auto
    (  long "num-leaves"
    <> metavar "NUMBER"
    <> help "Number of leaves in a tree"
    <> value 10
    <> showDefault
    )

data TraceState method_state = TraceState
  { trace_method_state :: !method_state
  }

trace
  :: Int
    -- ^ number of optimization steps
  -> Word64
    -- ^ random seed
  -> Int
    -- ^ number of observed sites
  -> Int
    -- ^ number of leaves in a tree
  -> IO ()
trace num_steps seed num_sites num_leaves =
  withFile "trace.csv" WriteMode $ \trace_h ->
  withFile "info.csv" WriteMode $ \info_h -> do
    hSetBuffering trace_h LineBuffering
    hSetBuffering info_h LineBuffering
    rnd_src <- newIORef (pureMT seed)

    (prob, true_bls) <- flip runRVar rnd_src $ do
      rate_mx <- rateMatrix bd
      topo <- topology bd
      bls <- branchLengths bd topo
      obs <- realisticObservations bd rate_mx topo bls
      return (Problem rate_mx obs topo, bls)

    let
      init_bls = 0.1 <$ true_bls
      true_ll = logLikelihood prob true_bls

    hPrintf info_h "var,value\n"
    hPrintf info_h "true_ll,%.6f\n" true_ll

    hPrintf trace_h "method,step,rate,ll,mse,time\n"

    forM_ methods $ \Method{..} -> do
      hPrintf trace_h "%s,0,NA,%.6f,%.6g,0\n"
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
              liftIO $ hPrintf trace_h "%s,%d,%.4g,%.6f,%.6g,%.3g\n" methodName istep actual_learning_rate ll (calculateMSE bls true_bls) time
              --liftIO $ print bls
            Nothing -> mzero

    return ()
  where
    bd = dnaBaseDistributions
      { numberOfSitesDistribution = pure num_sites
      , numberOfLeavesInTreeDistribution = pure num_leaves
      }
    methods =
      [noDescent] ++ (gradientDescent <$> [1e-4])

calculateMSE
  :: BranchLengths
  -> BranchLengths
  -> Double
calculateMSE bls true_bls =
  case sum . fmap (^(2::Int)) $ bls - true_bls of
    BranchLength sum_errors_squared ->
      sum_errors_squared / fromIntegral (size true_bls)
