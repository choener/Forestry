
{-# Options_GHC -fno-cse #-}

module Main where

import           Control.Lens
import           Data.ByteString.Lens
import           Control.Monad (void, when)
import           Data.Function ((&))
import qualified Data.Attoparsec.ByteString.Streaming as A
import qualified Data.ByteString.Streaming.Char8 as Q
import qualified Streaming as S
import qualified Streaming.Prelude as SP
import           System.Console.CmdArgs
import           System.Exit (exitSuccess, exitFailure)
import qualified Data.Vector as V
import           Numeric.Log

import           Biobase.Newick.Types
import           Biobase.Secondary.New (parseVienna, toTree, _label)
import           Biobase.Types.Structure
import           BioInf.ViennaRNA.RNAfold as RNAfold
import           Data.Forest.Static
import qualified Diagrams.TwoD.ProbabilityGrid as PG

import qualified Data.Forest.Static.Align.Affine as AA
import qualified Data.Forest.Static.Align.Linear as AL



-- | Which options do we provide.

data Options
  = AlignLinear
    { smatch        ∷ Int
    , smismatch     ∷ Int
    , sindel        ∷ Int
    , probabilities ∷ Bool
    , fillweight    ∷ PG.FillWeight
    , temperature   ∷ Double
    }
  | AlignAffine
    { smatch      ∷ Int
    , smismatch   ∷ Int
    , sindelopen  ∷ Int
    , sindelcont  ∷ Int
    , probabilities ∷ Bool
    , fillweight    ∷ PG.FillWeight
    , temperature   ∷ Double
    }
    deriving (Show,Data,Typeable)

oAlignLinear = AlignLinear
  { smatch        = 2         &= help "match score (2)"
  , smismatch     = (-2)      &= help "mismatch score (-2)"
  , sindel        = (-1)      &= help "in/del score (-1)"
  , probabilities = False     &= help "generate probability plot (it is advisable to pipe in RNAfold input with sequence identifiers)"
  , fillweight    = PG.FWlog  &= help ""
  , temperature   = 0.1
  }

oAlignAffine = AlignAffine
  { smatch = 2
  , smismatch = (-2)
  , sindelopen = (-1)
  , sindelcont = (-1)
  , probabilities = False
  , fillweight = PG.FWlog
  , temperature   = 0.1
  }

main ∷ IO ()
main = do
  o ← cmdArgs $ modes [oAlignLinear &= auto, oAlignAffine]
  case o of
    AlignLinear{} → runAlignLinear o
    AlignAffine{} → runAlignAffine o

-- | Read in pairs of @RNAfold@ outputs which are then fed to the linear
-- version of the tree alignment algorithm.

runAlignLinear ∷ Options → IO ()
runAlignLinear AlignLinear{..} = do
  let go ([]            S.:> r) = return r
      go ([(i,x)]       S.:> r) = return r
      go ([(i,x),(j,y)] S.:> r) = do
        let t1 = mkForest x
        let t2 = mkForest y
        AL.runAlignScoreTrees t1 t2 smatch smismatch sindel
        when probabilities $ do
          let fp' = x^.sequenceID.unpackedChars ++ y^.sequenceID.unpackedChars
          let fp = (++ "-" ++ show i ++ "-" ++ show j ++".eps") $ if null fp' then "prob" else fp'
          let t x = Exp $ fromIntegral x / temperature
          AL.runAlignScoreTreesIO
            fillweight PG.EPS
            fp
            t1 t2
            (t smatch) (t smismatch) (t sindel)
        return r
      go ∷ SP.Of [(Int,RNAfold)] a → IO a
  rest ← Q.getContents
       & A.parsed (RNAfold.pRNAfold RNAfold.ForceRNA 37)
       & void
       & SP.zip (SP.each [1..])
       & S.chunksOf 2
       & S.mapped SP.toList
       & S.mapsM_ go
  return ()

runAlignAffine ∷ Options → IO ()
runAlignAffine AlignAffine{..} = do
  let go ([]            S.:> r) = return r
      go ([(i,x)]       S.:> r) = return r
      go ([(i,x),(j,y)] S.:> r) = do
        let t1 = mkForest x
        let t2 = mkForest y
        AA.runAlignScoreTrees t1 t2 smatch smismatch sindelopen sindelcont
        when probabilities $ do
          let fp' = x^.sequenceID.unpackedChars ++ y^.sequenceID.unpackedChars
          let fp = (++ "-" ++ show i ++ "-" ++ show j ++".eps") $ if null fp' then "prob" else fp'
          let t x = Exp $ fromIntegral x / temperature
          AA.runAlignScoreTreesIO
            fillweight PG.EPS
            fp
            t1 t2
            (t smatch) (t smismatch) (t sindelopen) (t sindelcont)
        return r
      go ∷ SP.Of [(Int,RNAfold)] a → IO a
  rest ← Q.getContents
       & A.parsed (RNAfold.pRNAfold RNAfold.ForceRNA 37)
       & void
       & SP.zip (SP.each [1..])
       & S.chunksOf 2
       & S.mapped SP.toList
       & S.mapsM_ go
  return ()

mkForest ∷ RNAfold → Forest 'Pre V.Vector Info
mkForest x = forestPre [toTree (const $ Just $ Info "" 0) (Info "" 0) . either (error . show) id . parseVienna $ x^.mfe.foldedStructure.rnass]

