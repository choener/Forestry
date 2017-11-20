
module Data.Forest.Static.Align.Linear where

import           Control.Monad(forM_, unless)
import           Numeric.Log
import qualified Data.Text as Text
import qualified Data.Tree as T
import qualified Data.Vector as V
import qualified Data.Vector.Fusion.Stream.Monadic as SM
import qualified Data.Vector.Generic as VG
import           Text.Printf
import           Unsafe.Coerce

import           ADP.Fusion.Core
import           Biobase.Newick
import           Data.Forest.Static (TreeOrder(..),Forest)
import           Data.PrimitiveArray as PA hiding (map)
import           FormalLanguage.CFG
import qualified Data.Forest.Static as F
import qualified Diagrams.TwoD.ProbabilityGrid as PG

import           ADP.Fusion.Forest.Align.RL



[formalLanguage|
Verbose

Grammar: Global
N: T
N: F
N: M
T: n
S: [F,F]
[F,F] -> iter  <<< [T,T] [F,F]
[F,F] -> iter  <<< [M,M] [F,F]
[T,T] -> indel <<< [-,n] [F,F]
[T,T] -> delin <<< [n,-] [F,F]
[M,M] -> align <<< [n,n] [F,F]
[F,F] -> done  <<< [e,e]
//
Outside: Labolg
Source: Global
//

Emit: Global
Emit: Labolg
|]

makeAlgebraProduct ''SigGlobal

resig :: Monad m => SigGlobal m a b c d -> SigLabolg m a b c d
resig (SigGlobal gdo git gal gin gde gh) = SigLabolg gdo git gal gin gde gh
{-# Inline resig #-}

score :: Monad m => Int -> Int -> Int -> SigGlobal m Int Int Info Info
score matchSc notmatchSc delinSc = SigGlobal
  { gDone  = \ (Z:.():.())  -> 0
  , gIter  = \ t f          -> t+f
  , gAlign = \ (Z:.a:.b) f  -> f + if label a == label b then matchSc else notmatchSc
  , gIndel = \ (Z:.():.b) f -> f + delinSc
  , gDelin = \ (Z:.a:.()) f -> f + delinSc
  , gH     = SM.foldl' max (-88888)
  }
{-# Inline score #-}

part :: Monad m => Log Double -> Log Double -> Log Double -> Log Double -> SigGlobal m (Log Double) (Log Double) Info Info
part matchSc mismatchSc indelSc temp = SigGlobal
  { gDone  = \ (Z:.():.())  -> 1
  , gIter  = \ t f          -> t * f
  , gAlign = \ (Z:.a:.b) f  -> f * if label a == label b then matchSc else mismatchSc
  , gIndel = \ (Z:.():.b) f -> f * indelSc
  , gDelin = \ (Z:.a:.()) f -> f * indelSc
  , gH     = SM.foldl' (+) 0.0000001
  }
{-# Inline part #-}

type Pretty = [[T.Tree (Info,Info)]]
pretty :: Monad m => SigGlobal m [T.Tree (Info,Info)] [[T.Tree ((Info,Info))]] Info Info
pretty = SigGlobal
  { gDone  = \ (Z:.():.())  -> []
  , gIter  = \ t f          -> t++f
  , gAlign = \ (Z:.a:.b) f  -> [T.Node (a,b) f]
  , gIndel = \ (Z:.():.b) f -> [T.Node (Info "-" 0,b) f]
  , gDelin = \ (Z:.a:.()) f -> [T.Node (a,Info "-" 0) f]
  , gH     = SM.toList
  }
{-# Inline pretty #-}



type Trix = TreeIxR Pre V.Vector Info I
type Tbl x = TwITbl Id Unboxed (Z:.EmptyOk:.EmptyOk) (Z:.Trix:.Trix) x
type Frst = Forest Pre V.Vector Info
type TblBt x = TwITblBt Unboxed (Z:.EmptyOk:.EmptyOk) (Z:.Trix:.Trix) Int Id Id [x]
type B = T.Tree (Info,Info)

runForward :: Frst -> Frst -> Int -> Int -> Int -> Z:.Tbl Int :.Tbl Int:.Tbl Int
runForward f1 f2 matchSc notmatchSc delinSc = mutateTablesDefault $
                   gGlobal (score matchSc notmatchSc delinSc)
                   (ITbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99999) [] ))
                   (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99998) [] ))
                   (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99997) [] ))
                   (node NTany $ F.label f1)
                   (node NTany $ F.label f2)
{-# NoInline runForward #-}

runInside :: Frst -> Frst -> Log Double -> Log Double -> Log Double -> Log Double -> Z:.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double)
runInside f1 f2 matchSc mismatchSc indelSc temperature = mutateTablesDefault $
                   gGlobal (part matchSc mismatchSc indelSc temperature)
                   (ITbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000001) [] ))
                   (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000002) [] ))
                   (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000003) [] ))
                   (node NTany $ F.label f1)
                   (node NTany $ F.label f2)
{-# NoInline runInside #-}

type Trox = TreeIxR Pre V.Vector Info O
type OTbl x = TwITbl Id Unboxed (Z:.EmptyOk:.EmptyOk) (Z:.Trox:.Trox) x

runOutside :: Frst -> Frst -> Log Double -> Log Double -> Log Double -> Log Double -> Z:.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double) -> Z:.OTbl (Log Double):.OTbl (Log Double):.OTbl (Log Double)
runOutside f1 f2 matchSc mismatchSc indelSc temperature (Z:.iF:.iM:.iT)
  = mutateTablesDefault $
    gLabolg (resig (part matchSc mismatchSc indelSc temperature))
    (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000001) [] ))
    (ITbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000002) [] ))
    (ITbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000003) [] ))
    iF
    iM
    iT
    (node NTany $ F.label f1)
    (node NTany $ F.label f2)
{-# NoInline runOutside #-}

runS :: Frst -> Frst -> Int -> Int -> Int -> (Z:.Tbl Int :.Tbl Int:.Tbl Int, Int ,[[T.Tree (Info, Info)]] )
runS f1 f2 matchSc notmatchSc delinSc = (fwd,unId $ axiom f, unId $ axiom fb)
  where fwd@(Z:.f:.m:.t) = runForward f1 f2 matchSc notmatchSc delinSc
        Z:.fb:.fm:.tb = gGlobal ((score matchSc notmatchSc delinSc) <|| pretty) (toBacktrack f (undefined :: Id a -> Id a)) (toBacktrack m (undefined :: Id a -> Id a)) (toBacktrack t (undefined :: Id a -> Id a))
                        (node NTany $ F.label f1) (node NTany $ F.label f2)
                        :: Z:.TblBt B:.TblBt B:.TblBt B
{-# NoInline runS #-}

runIO f1 f2 matchSc mismatchSc indelSc temperature = (fwd,out,unId $ axiom f)
  where fwd@(Z:.f:.m:.t) = runInside f1 f2 matchSc mismatchSc indelSc temperature
        out@(Z:.oft:.omt:.ott) = runOutside f1 f2 matchSc mismatchSc indelSc temperature fwd
{-# NoInline runIO #-}

--         a            a
--        / \          / \
--       e   d        b   f
--      / \              / \
--     b   c            c   d
--
--                  (a,a)                 100
--              /          \
--         (e,-)            (-,f)         (-3) (-5)
--        /     \          /     \
--   (b,b)       (c,-) (-,c)      (d,d)   100  (-3) (-5) 100
--
--
--
--             (a,a)                          100
--            /     \
--       (e,-)       (d,-)                    (-3) (-3)
--      /     \
-- (b,b)       (-,f)                          100  (-5)
--            /     \
--       (c,c)       (-,d)                    100  (-5)

t11 = "a;"
t12 = "a;"
t21 = "(b,c)a;"
t22 = "(b,c)a;"
t31 = "((d,e,f)b,(z)c)a;"  --
t32 = "(((d,e)y,f)b,(c,(x)i)g)a;"  --
t41 = "d;(b)e;" -- (b,c)e;"    -- '-3'
t42 = "(d)f;b;" -- b;"
t51 = "(b:1,c:1)a:1;"
t52 = "b:2;c:2;"
t61 = "((b,c)e,d)a;"
t62 = "(b,(c,d)f)a;"
t71 = "(b)a;"
t72 = "(b)a;"


test f1 f2 = do
  n1 ← readFile f1
  n2 ← readFile f2
  runAlignS n1 n2 1 (-1) (-1)
{-# NoInline test #-}

runAlignScoreTrees t1 t2 matchSc notmatchSc delinSc = do
  let (fwd,sc,bt') = runS t1 t2 matchSc notmatchSc delinSc
  let (Z:.TW (ITbl _ _ _ ift) _ :. TW (ITbl _ _ _ imt) _ :. TW (ITbl _ _ _ itt) _) = fwd
  let bt = take 1 bt' -- TODO make nice !!! nub bt'
  printf "Score: %10d\n" sc
  forM_ bt $ \b -> do
    putStrLn ""
    forM_ b $ \x -> putStrLn $ T.drawTree $ fmap show x
{-# NoInline runAlignScoreTrees #-}

runAlignS t1' t2' matchSc notmatchSc delinSc = do
  let f x = either error (F.forestPre . map getNewickTree) $ newicksFromText x
      t1 = f $ Text.pack t1'
      t2 = f $ Text.pack t2'
  let (fwd,sc,bt') = runS t1 t2 matchSc notmatchSc delinSc
  let (Z:.TW (ITbl _ _ _ ift) _ :. TW (ITbl _ _ _ imt) _ :. TW (ITbl _ _ _ itt) _) = fwd
  let bt = take 1 bt' -- TODO make nice !!! nub bt'
  printf "Score: %10d\n" sc
  forM_ bt $ \b -> do
    putStrLn ""
    forM_ b $ \x -> putStrLn $ T.drawTree $ fmap show x
{-# NoInline runAlignS #-}

runAlignIO fw probFileTy probFile t1' t2' matchSc mismatchSc indelSc temperature = do
  let f x = either error (F.forestPre . map getNewickTree) $ newicksFromText x
      t1 = f $ Text.pack t1'
      t2 = f $ Text.pack t2'
  let (inn,out,_) = runIO t1 t2 matchSc mismatchSc indelSc temperature -- (t2 {F.lsib = VG.fromList [-1,-1], F.rsib = VG.fromList [-1,-1]})
  let (Z:.TW (ITbl _ _ _ ift) _ :. TW (ITbl _ _ _ imt) _ :. TW (ITbl _ _ _ itt) _) = inn
  let (Z:.TW (ITbl _ _ _ oft) _ :. TW (ITbl _ _ _ omt) _ :. TW (ITbl _ _ _ ott) _) = out
  let (Z:.(TreeIxR frst1 lb1 _):.(TreeIxR frst2 lb2 _), Z:.(TreeIxR _ ub1 _):.(TreeIxR _ ub2 _)) = bounds oft
  let ix = (Z:.TreeIxR frst1 lb1 F:.TreeIxR frst2 lb2 F)
  let scift = ift ! ix
  print scift
  let scoft = Prelude.sum [ oft ! (Z:.TreeIxR frst1 b1 F :. TreeIxR frst2 b2 F) | b1 <- [lb1 .. ub1], b2 <- [lb2 .. ub2] ]
  print scoft
  let scimt = Prelude.sum [ imt ! (Z:.TreeIxR frst1 b1 T :. TreeIxR frst2 b2 T) | b1 <- [lb1 .. ub1], b2 <- [lb2 .. ub2] ]
  print scimt
  let scomt = Prelude.sum [ omt ! (Z:.TreeIxR frst1 b1 T :. TreeIxR frst2 b2 T) | b1 <- [lb1 .. ub1], b2 <- [lb2 .. ub2] ]
  print scomt
  let ps = map (\(k,k1,k2) ->
            let k' = unsafeCoerce k
            in  ( k1
                , k2
                , ((imt!k) * (omt!k') / scift)
                , (maybe "-" label $ F.label t1 VG.!? k1)
                , (maybe "-" label $ F.label t2 VG.!? k2)
                )) [ (Z:.TreeIxR frst1 k1 T:.TreeIxR frst2 k2 T,k1,k2) | k1 <- [lb1 .. ub1 - 1], k2 <- [lb2 .. ub2 - 1] ]
  --
  let gsc = map (\(k1,k2,sc,l1,l2) -> sc) ps
  let fillText [] = " "
      fillText xs = xs
  let gl1 = map (\k1 -> fillText . Text.unpack $ (maybe "-" label $ F.label t1 VG.!? k1)) [lb1 .. ub1 - 1]
  let gl2 = map (\k2 -> fillText . Text.unpack $ (maybe "-" label $ F.label t2 VG.!? k2)) [lb2 .. ub2 - 1]
  case probFileTy of
         PG.SVG -> PG.svgGridFile probFile fw PG.FSfull ub1 ub2 gl1 gl2 gsc
         PG.EPS -> PG.epsGridFile probFile fw PG.FSfull ub1 ub2 gl1 gl2 gsc
{-# NoInline runAlignIO #-}

