{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE RecordWildCards #-}

-- | Compressive Sampling Matched Pursuit (CoSaMP): Iterative signal recovery
-- from incomplete and inaccurate samples" by Deanna Needell & Joel Tropp
module CoSaMP ( cosamp
              , dctOrtho
              , idctOrtho
              ) where

import Numeric.LinearAlgebra
import Numeric.LinearAlgebra.Devel        (ti)
import Prelude                     hiding ((<>), reverse)
import Data.List                          (nub)
import Control.Monad.State
import Statistics.Transform               (dct, idct)

-- | CoSaMP Algorithm State
data CoSaMPState = CoSaMPState 
    { p :: Double        -- ^ Precision
    , i :: Int           -- ^ Iteration
    , e :: Double        -- ^ Tolerance
    , φ :: Matrix Double -- ^ Sampling Matrix
    , u :: Vector Double -- ^ Noisy sampling vector
    , a :: Vector Double -- ^ s-sparse approximation of target signal @ iteration @i@
    , v :: Vector Double -- ^ Updated sample vector
    , s :: Int           -- ^ Sparsity level
    } deriving (Show)

-- | Compressive Sampling Matched Pursuit (CoSaMP)
cosamp :: Matrix Double -- ^ sampling matrix Φ
       -> Vector Double -- ^ noisy sample vector u
       -> Int           -- ^ sparsity level s
       -> Int           -- ^ Number of Iterations (default 1000)
       -> Double        -- ^ Tolerance (default 1.0e-10)
       -> Vector Double -- ^ s-sparse approsimation of the target signal
cosamp φ' u' s' i' e' = evalState (cosamp' False) st
  where
    p' = 1.0e-12
    a' = konst (0 :: Double) . snd $ size φ'
    st = CoSaMPState { p = p'
                     , i = i'
                     , e = e'
                     , φ = φ'
                     , u = u'
                     , a = a'
                     , v = u'
                     , s = s' }

-- | CoSaMP Algorithm in the CoSaMPState Monad
cosamp' :: Bool -> State CoSaMPState (Vector Double)
cosamp' True  = gets a
cosamp' False = do
    st@CoSaMPState{..} <- get
    let y  = abs $ tr φ #> v
        c  = sortVector y ! (size y - 1 - 2*s)
        ω  = find (\val -> val > c && val > p) y
        t  = nub $ ω ++ find (/= 0) a :: [IndexOf Vector]
        b  = pinv (φ ?? (All, Pos (idxs t))) #> u
        g  = sortVector (abs b) ! (size b - 1 - s)
        j  = idxs . find (\b'' -> b'' > g && b'' > p) $ abs b
        t' = flatten $ asRow (idxs t) ?? (All, Pos j)
        b' = flatten $ asRow b ?? (All, Pos j)
        a' = accum a const $ zip (map ti $ toList t') (toList b')
        v' = u - (φ ?? (All, Pos t') #> b')
        i' = i - 1
        h  = (norm_2 v / norm_2 u) < e || i <= 0
    put $ st { i = i', v = v', a = a' }
    cosamp' h

-- | DCT Type II where corresponding matrix coefficients are made orthonormal
dctOrtho :: Matrix Double -> Matrix Double
dctOrtho x = y 
  where
    n  = fst $ size x
    n' = realToFrac n :: Double
    t  = fromColumns . map dct $ toColumns x
    r0 = scale (1.0 / ( 2.0 * sqrt n')) (t ?? (Pos $ idxs [0], All))
    rs = scale (sqrt (2.0 / n') / 2.0) (t ?? (Pos $ idxs [1 .. (n - 1)], All))
    y  = r0 === rs

-- | Inverse DCT Type II where corresponding matrix coefficients are made orthonormal
idctOrtho :: Vector Double -> Vector Double
idctOrtho y = x
  where
    n = realToFrac $ size y :: Double
    x = scale (0.5 * sqrt (2.0 / n)) (idct y) 
