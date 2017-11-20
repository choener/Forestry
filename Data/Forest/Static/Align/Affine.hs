
module Data.Forest.Static.Align.Affine where

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
N: T -- tree
N: F -- forest
N: Z -- tree for gaps
N: Q -- sibling gap mode
N: R -- parent gap mode
N: E
T: n
S: [F,F]
[F,F] -> iter    <<< [T,T] [F,F]
[F,F] -> fgap    <<< [T,Z] [Q,Q]
[F,F] -> fgap    <<< [Z,T] [Q,Q]
[Z,T] -> indel   <<< [-,n] [R,R]
[T,Z] -> delin   <<< [n,-] [R,R]
[T,T] -> align   <<< [n,n] [F,F]
[F,F] -> done    <<< [E,E]
[R,R] -> done    <<< [E,E]
[R,R] -> pgap <<< [T,T] [R,R]
[R,R] -> pgap <<< [T,Z] [R,R]
[R,R] -> pgap <<< [Z,T] [R,R]
[Q,Q] -> done    <<< [E,E]
[Q,Q] -> siter <<< [T,T] [F,F]
[Q,Q] -> sgap <<< [T,Z] [Q,Q]
[Q,Q] -> sgap <<< [Z,T] [Q,Q]
[E,E] -> finalDone <<< [e,e]
//
Outside: Labolg
Source: Global
//

Emit: Global
Emit: Labolg
|]

makeAlgebraProduct ''SigGlobal

resig :: Monad m => SigGlobal m a b c d -> SigLabolg m a b c d
resig (SigGlobal gdo git gsi gal gin gde gfg gpg gsg gfi gh) = SigLabolg gdo git gsi gal gin gde gfg gpg gsg gfi gh
{-# Inline resig #-}


score :: Monad m => Int -> Int -> Int -> Int -> SigGlobal m Int Int Info Info
score matchSc notmatchSc delinSc affinSc = SigGlobal -- match affine deletion 
  { gDone  = \ f -> f --  (Z:.():.()) -> 0 -- traceShow "EEEEEEEEEEEEE" 0
  , gFinalDone  = \ (Z:.():.()) -> 0 -- traceShow "EEEEEEEEEEEEE" 0
  , gIter  = \ t f -> t+f
  , gSiter  = \ t f -> t+f
  , gAlign = \ (Z:.c:.b) f -> tSI glb ("ALIGN",f,c,b) $ f + if label c == label b then matchSc else notmatchSc
  , gIndel = \ (Z:.():.b) f -> tSI glb ("INDEL",f,b) $ f
  , gDelin = \ (Z:.c:.()) f -> tSI glb ("DELIN",f,c) $ f
  , gFgap = \ t f -> tSI glb ("gap",f+t,delinSc) $ t + f + delinSc --gap open
  , gPgap = \ t f -> tSI glb ("gap",f+t,affinSc) $ t + f + affinSc --gap extension
  , gSgap = \ t f -> tSI glb ("gap",f+t,affinSc) $ t + f + affinSc --gap extension
  , gH     = SM.foldl' max (-88888)
  }
{-# Inline score #-}

part :: Monad m => Log Double -> Log Double -> Log Double -> Log Double -> Log Double -> SigGlobal m (Log Double) (Log Double) Info Info
part matchSc mismatchSc delinSc affinSc temp = SigGlobal
  { gDone  = \ f -> f 
  , gIter  = \ t f -> tSI glb ("TFTFTFTFTF",t,f) $ t * f
  , gFinalDone  = \ (Z:.():.()) -> 1
  , gSiter  = \ t f -> tSI glb ("TFTFTFTFTF",t,f) $ t * f
  , gAlign = \ (Z:.a:.b) f -> tSI glb ("ALIGN",f,a,b) $ f * if label a == label b then matchSc else mismatchSc
  , gIndel = \ (Z:.():.b) f -> tSI glb ("INDEL",f,b) $ f
  , gDelin = \ (Z:.a:.()) f -> tSI glb ("DELIN",f,a) $ f
  , gFgap = \ t f -> t * f * delinSc
  , gPgap = \ t f -> t * f * affinSc
  , gSgap = \ t f -> t * f * affinSc
  , gH     = SM.foldl' (+) 0.000000
  }
{-# Inline part #-}



type Pretty' = [[T.Tree (Info,Info)]]
pretty' :: Monad m => SigGlobal m [T.Tree (Info,Info)] [[T.Tree ((Info,Info))]] Info Info
pretty' = SigGlobal
  { gDone  = \ f -> f -- (Z:.():.()) -> []
  , gFinalDone  = \ (Z:.():.()) -> []
  , gIter  = \ t f -> t++f
  , gSiter  = \ t f -> t++f
  , gAlign = \ (Z:.a:.b) f -> [T.Node (a,b) f]
  , gIndel = \ (Z:.():.b) f -> [T.Node (Info "-" 0,b) f]
  , gDelin = \ (Z:.a:.()) f -> [T.Node (a,Info "-" 0) f]
  , gPgap = \ t f -> t ++ f
  , gSgap = \ t f -> t ++ f
  , gFgap = \ t f -> t ++ f
  , gH     = SM.toList
  }
{-# Inline pretty' #-}



type Trix = TreeIxR Pre V.Vector Info I
type Tbl x = TwITbl Id Unboxed (Z:.EmptyOk:.EmptyOk) (Z:.Trix:.Trix) x
type Frst = Forest Pre V.Vector Info
type TblBt x = TwITblBt Unboxed (Z:.EmptyOk:.EmptyOk) (Z:.Trix:.Trix) Int Id Id [x]
type B = T.Tree (Info,Info)


-- | likelihood part
--
--NOTE for an explanation which ITbl gets @0@ or @1@ check @runInside@,
--they have the same order requirements.

runForward :: Frst -> Frst -> Int -> Int -> Int -> Int -> Z:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int
runForward f1 f2 matchSc notmatchSc delinSc affinSc = let
                         in
                           mutateTablesDefault $
                           gGlobal (score matchSc notmatchSc delinSc affinSc) -- costs
                           (iTbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99999) [] ))
                           (iTbl 0 2 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99998) [] ))
                           (iTbl 0 2 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99997) [] ))
                           (iTbl 0 2 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99996) [] ))
                           (iTbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99995) [] ))
                           (iTbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99994) [] ))
                           (iTbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (-99993) [] ))
                           (node NTany $ F.label f1)
                           (node NTany $ F.label f2)
{-# NoInline runForward #-}

iTbl = ITbl
{-# NoInline iTbl #-}

-- |inside part
--
-- NOTE : Each @ITbl@ has a big-order (1st index) and a little-order index.
-- We need these indices because we operate on unboxed tables internally
-- for efficiency reasons.
--
-- The big-order index is for full table calculations. All @ITbl@s with
-- smaller "big order" are fully calculated before any @ITbl@s with larger
-- big orders.
--
-- The little-order influences per-cell calculations. Cells with the same
-- index are sorted by order and those with a smaller index calculated
-- first.
--
-- For this calculation we have the following:
--
-- @EE@ has only one rule: @EE -> ε@, the full table is filled before any
-- other table and has big order @0@.
--
-- Now, we only have tables with big order @1@ left.
--
-- @TT@, @TZ@, and @ZT@ have only rules that each have at least one
-- non-empty terminal in play. @TZ@ for example has @TZ -> [n,-] RR@.
-- They are executed first for a given index, because they are guaranteed
-- to get "smaller" during recursion. Hence little order @0@.
--
-- Finally, @FF@, @QQ@, @RR@ all combine *only* syntactic symbols on the
-- right-hand side. In addition these symbols can be empty. They need to
-- come last and have little order @1@.

runInside :: Frst -> Frst -> Log Double -> Log Double -> Log Double -> Log Double -> Log Double -> Z:.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double)
runInside f1 f2 matchSc notmatchSc delinSc affinSc temperature = let
                         in
                           mutateTablesDefault $
                           gGlobal (part matchSc notmatchSc delinSc affinSc temperature) -- costs
                           (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000001) [] ))   -- EE
                           (ITbl 1 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000002) [] ))   -- FF
                           (ITbl 1 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000003) [] ))   -- QQ
                           (ITbl 1 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000004) [] ))   -- RR
                           (ITbl 1 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000005) [] ))   -- TT
                           (ITbl 1 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000006) [] ))   -- TZ
                           (ITbl 1 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000007) [] ))   -- ZT
                           (node NTany $ F.label f1)
                           (node NTany $ F.label f2)
{-# NoInline runInside #-}



-- outside part
type Trox = TreeIxR Pre V.Vector Info O
type OTbl x = TwITbl Id Unboxed (Z:.EmptyOk:.EmptyOk) (Z:.Trox:.Trox) x

-- | Actual outside calculations.
--
-- NOTE: The big and little order indices are reversed compared to
-- @runInside@. We need to evaluate the tables outside in now.

runOutside :: Frst -> Frst -> Log Double -> Log Double -> Log Double -> Log Double -> Log Double -> Z:.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double):.Tbl (Log Double) -> Z:.OTbl (Log Double):.OTbl (Log Double):.OTbl (Log Double):.OTbl (Log Double):.OTbl (Log Double):.OTbl (Log Double):.OTbl (Log Double)
runOutside f1 f2 matchSc mismatchSc indelSc affinSc temperature (Z:.iE:.iF:.iQ:.iR:.iT:.iS:.iZ)
  = mutateTablesDefault $
    gLabolg (resig (part matchSc mismatchSc indelSc affinSc temperature))
    (ITbl 1 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000001) [] ))
    (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000002) [] ))
    (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000003) [] ))
    (ITbl 0 0 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000004) [] ))
    (ITbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000005) [] ))
    (ITbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000006) [] ))
    (ITbl 0 1 (Z:.EmptyOk:.EmptyOk) (PA.fromAssocs (Z:.minIx f1:.minIx f2) (Z:.maxIx f1:.maxIx f2) (0.00000007) [] ))
    iF
    iQ
    iR
    iT
    iS
    iZ
    (node NTany $ F.label f1)
    (node NTany $ F.label f2)
{-# NoInline runOutside #-}



-- inside part
run :: Frst -> Frst -> Int -> Int -> Int -> Int -> (Z:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int:.Tbl Int,Int,Pretty')
run f1 f2 matchSc notmatchSc delinSc affinSc = (fwd,unId $ axiom a2, unId $ axiom b2)
  where fwd@(Z:.a1:.a2:.a3:.a4:.a5:.a6:.a7) = runForward f1 f2 matchSc notmatchSc delinSc affinSc
        Z:.b1:.b2:.b3:.b4:.b5:.b6:.b7 
                    = gGlobal ((score matchSc notmatchSc delinSc affinSc) <|| pretty') 
                    (toBacktrack a1 (undefined :: Id a -> Id a)) 
                    (toBacktrack a2 (undefined :: Id a -> Id a))  
                    (toBacktrack a3 (undefined :: Id a -> Id a))  
                    (toBacktrack a4 (undefined :: Id a -> Id a))  
                    (toBacktrack a5 (undefined :: Id a -> Id a))  
                    (toBacktrack a6 (undefined :: Id a -> Id a))  
                    (toBacktrack a7 (undefined :: Id a -> Id a))  
                    (node NTany $ F.label f1) (node NTany $ F.label f2)
                    :: Z:.TblBt B:.TblBt B:.TblBt B:.TblBt B:.TblBt B:.TblBt B:.TblBt B
{-# NoInline run #-}



-- outside part

runIO f1 f2 matchSc mismatchSc indelSc affinSc temperature = (fwd,out,unId $ axiom f)
  where fwd@(Z:.e:.f:.q:.r:.t:.s:.z) = runInside f1 f2 matchSc mismatchSc indelSc affinSc temperature
        out@(Z:.oet:.oft:.oqt:.ort:.ott:.ost:.ozt) = runOutside f1 f2 matchSc mismatchSc indelSc affinSc temperature fwd
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

-- all new from here


test n1 n2 = do
  f1 <- readFile n1
  f2 <- readFile n2
  runAlignS f1 f2 10 (-30) (-10) (-1)
{-# NoInline test #-}

runAlignScoreTrees t1 t2 matchSc notmatchSc delinSc affinSc = do
  let (fwd,sc,bt') = run t1 t2 matchSc notmatchSc delinSc affinSc
  let (Z:.TW (ITbl _ _ _ iet) _ :.TW (ITbl _ _ _ ift) _ :.TW (ITbl _ _ _ iqt) _ :.TW (ITbl _ _ _ irt) _ :.TW (ITbl _ _ _ itt) _ :.TW (ITbl _ _ _ ist) _ :.TW (ITbl _ _ _ izt) _) = fwd
  let bt = take 1 bt' -- nub bt'
  printf "Score: %10d\n" sc
  forM_ bt $ \b -> do
    putStrLn ""
    forM_ b $ \x -> putStrLn $ T.drawTree $ fmap show x
{-# NoInline runAlignScoreTrees #-}

runAlignS t1' t2' matchSc notmatchSc delinSc affinSc = do
  let f x = either error (F.forestPre . map getNewickTree) $ newicksFromText x
      t1 = f $ Text.pack t1'
      t2 = f $ Text.pack t2'
  let (fwd,sc,bt') = run t1 t2 matchSc notmatchSc delinSc affinSc
  let (Z:.TW (ITbl _ _ _ iet) _ :.TW (ITbl _ _ _ ift) _ :.TW (ITbl _ _ _ iqt) _ :.TW (ITbl _ _ _ irt) _ :.TW (ITbl _ _ _ itt) _ :.TW (ITbl _ _ _ ist) _ :.TW (ITbl _ _ _ izt) _) = fwd
  let bt = take 1 bt' -- nub bt'
  printf "Score: %10d\n" sc
--  forM_ bt $ \b -> do
--    putStrLn ""
--    forM_ b $ \x -> putStrLn $ T.drawTree $ fmap show x
{-# NoInline runAlignS #-}

runAlignIO fw probFileTy probFile t1' t2' matchSc mismatchSc indelSc affinSc temperature = do
  let f x = either error (F.forestPre . map getNewickTree) $ newicksFromText x
      t1 = f $ Text.pack t1'
      t2 = f $ Text.pack t2'
  let (inn,out,_) = runIO t1 t2 matchSc mismatchSc indelSc affinSc temperature -- (t2 {F.lsib = VG.fromList [-1,-1], F.rsib = VG.fromList [-1,-1]})
  let (Z:.TW (ITbl _ _ _ iet) _ :.TW (ITbl _ _ _ ift) _ :.TW (ITbl _ _ _ iqt) _ :.TW (ITbl _ _ _ irt) _ :.TW (ITbl _ _ _ itt) _ :.TW (ITbl _ _ _ ist) _ :.TW (ITbl _ _ _ izt) _) = inn
  let (Z:.TW (ITbl _ _ _ iet) _ :.TW (ITbl _ _ _ oft) _ :.TW (ITbl _ _ _ oqt) _ :.TW (ITbl _ _ _ ort) _ :.TW (ITbl _ _ _ ott) _ :.TW (ITbl _ _ _ ost) _ :.TW (ITbl _ _ _ ozt) _) = out
  let (Z:.(TreeIxR frst1 lb1 _):.(TreeIxR frst2 lb2 _), Z:.(TreeIxR _ ub1 _):.(TreeIxR _ ub2 _)) = bounds oft
  let ix = (Z:.TreeIxR frst1 lb1 F:.TreeIxR frst2 lb2 F)
  let sc = ift ! ix
  print sc
  let ps = map (\(k,k1,k2) ->
            let k' = unsafeCoerce k
            in  ( k1
                , k2
                , ((itt!k) * (ott!k') / sc)
                , (maybe "-" label $ F.label t1 VG.!? k1)
                , (maybe "-" label $ F.label t2 VG.!? k2)
                )) [ (Z:.TreeIxR frst1 k1 T:.TreeIxR frst2 k2 T,k1,k2) | k1 <- [lb1 .. ub1 - 1], k2 <- [lb2 .. ub2 - 1] ]
  --
  let gsc = map (\(k1,k2,sc,l1,l2) -> sc) ps
  let fillText [] = " "
      fillText xs = xs
  mapM_ print gsc
  let gl1 = map (\k1 -> fillText . Text.unpack $ (maybe "-" label $ F.label t1 VG.!? k1)) [lb1 .. ub1 - 1]
  let gl2 = map (\k2 -> fillText . Text.unpack $ (maybe "-" label $ F.label t2 VG.!? k2)) [lb2 .. ub2 - 1]
  case probFileTy of
         PG.SVG -> PG.svgGridFile probFile fw PG.FSfull ub1 ub2 gl1 gl2 gsc
         PG.EPS -> PG.epsGridFile probFile fw PG.FSfull ub1 ub2 gl1 gl2 gsc
{-# NoInline runAlignIO #-}

