
{-# Options_GHC -fno-cse #-}

module Main where

import           Control.Lens
import           Control.Monad (void)
import           Data.Function ((&))
import qualified Data.Attoparsec.ByteString.Streaming as A
import qualified Data.ByteString.Streaming.Char8 as Q
import qualified Streaming as S
import qualified Streaming.Prelude as SP
import           System.Console.CmdArgs
import           System.Exit (exitSuccess, exitFailure)
import qualified Data.Vector as V

import           BioInf.ViennaRNA.RNAfold as RNAfold
import           Biobase.Secondary.New (parseVienna, toTree, _label)
import           Data.Forest.Static
import           Biobase.Types.Structure
import           Biobase.Newick.Types

import qualified Data.Forest.Static.Align.Affine as AA
import qualified Data.Forest.Static.Align.Linear as AL



-- | Which options do we provide.

data Options
  = AlignLinear
    { smatch        ∷ Int
    , smismatch     ∷ Int
    , sindel        ∷ Int
    , probabilities ∷ Bool
    }
  | AlignAffine
    { smatch      ∷ Int
    , smismatch   ∷ Int
    , sindelopen  ∷ Int
    , sindelcont  ∷ Int
    }
    deriving (Show,Data,Typeable)

oAlignLinear = AlignLinear
  { smatch        = 2     &= help "match score (2)"
  , smismatch     = (-2)  &= help "mismatch score (-2)"
  , sindel        = (-1)  &= help "in/del score (-1)"
  , probabilities = False &= help "generate probability plot (it is advisable to pipe in RNAfold input with sequence identifiers)"
  }

oAlignAffine = AlignAffine
  { smatch = 2
  , smismatch = (-2)
  , sindelopen = (-1)
  , sindelcont = (-1)
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
  let go ([]    S.:> r) = return r
      go ([x]   S.:> r) = return r
      go ([x,y] S.:> r) = do
        let t1 = mkForest x
        let t2 = mkForest y
        AL.runAlignScoreTrees t1 t2 smatch smismatch sindel
        return r
      go ∷ SP.Of [RNAfold] a → IO a
  rest ← Q.getContents
       & A.parsed (RNAfold.pRNAfold RNAfold.ForceRNA 37)
       & S.chunksOf 2
       & S.mapped SP.toList
       & S.mapsM_ go
  return ()
--  case rest of
--    Left (msg, bla) → do
--      print msg
--      exitFailure
--    Right () → exitSuccess

runAlignAffine ∷ Options → IO ()
runAlignAffine AlignAffine{..} = do
  let go ([]    S.:> r) = return r
      go ([x]   S.:> r) = return r
      go ([x,y] S.:> r) = do
        let t1 = mkForest x
        let t2 = mkForest y
        AA.runAlignScoreTrees t1 t2 smatch smismatch sindelopen sindelcont
        return r
      go ∷ SP.Of [RNAfold] a → IO a
  rest ← Q.getContents
       & A.parsed (RNAfold.pRNAfold RNAfold.ForceRNA 37)
       & S.chunksOf 2
       & S.mapped SP.toList
       & S.mapsM_ go
  return ()

mkForest ∷ RNAfold → Forest 'Pre V.Vector Info
mkForest x = forestPre [toTree (const $ Just $ Info "" 0) (Info "" 0) . either (error . show) id . parseVienna $ x^.mfe.foldedStructure.rnass]

