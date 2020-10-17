import Options.Applicative
import Control.Monad
import Control.Monad.State
import Data.Random
import Text.Printf (printf)
import Phylogenetics.Types
import Phylogenetics.Gen
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
  , trace_bls :: !BranchLengths
  }

trace :: Int -> IO ()
trace num_steps = do
  (prob, true_bls) <- sample $ do
    rate_mx <- rateMatrix bd
    topo <- topology bd
    bls <- branchLengths bd topo
    obs <- realisticObservations bd rate_mx topo bls
    return (Problem rate_mx obs topo, bls)

  let init_bls = 0.1 <$ true_bls

  printf "method,step,rate,ll,mse\n"

  forM_ methods $ \Method{..} ->
    flip evalStateT (TraceState methodInit init_bls) $
      forM_ [1 .. num_steps] $ \istep -> do
        prev_state <- get
        case methodStep prob (trace_bls prev_state) (trace_method_state prev_state) of
          Just (bls, actual_learning_rate, ll, st) -> do
            put $! TraceState st bls
            let mse =
                  case sum . fmap (^(2::Int)) $ bls - true_bls of
                    BranchLength sum_errors_squared ->
                      sum_errors_squared / fromIntegral (size true_bls)
            liftIO $ printf "%s,%d,%.4g,%.6f,%.6g\n" methodName istep actual_learning_rate ll mse
          Nothing -> mzero

  return ()
  where
    bd = dnaBaseDistributions
    methods = gradientDescent <$> [1, 0.1, 0.01]
