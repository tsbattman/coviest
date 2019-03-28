{-# LANGUAGE ViewPatterns #-}

module Statistics.Covariance.Shrinkage.MMSE (
    mmseCov
  , shrinkCov
  , sampleCov
  , naiveCov
  , pinCoef
  , lwShrinkageCoef
  , rblwShrinkCovCoef
  , rblwShrinkageCoef
  , oracleShrinkCovCoef
  , oracleShrinkageCoef
  ) where

import Numeric.LinearAlgebra hiding ((<>))
import Numeric.Sum (kbn, sumVector)
import Statistics.Sample (mean)
import qualified Data.Vector.Unboxed as VU

mmseCov :: (Matrix Double -> Double) -> Matrix Double -> Herm Double
mmseCov shrink x = shrinkCov s p f
  where
    s = sampleCov x
    f = naiveCov x
    p = shrink x

shrinkCov :: Herm Double -> Double -> Herm Double -> Herm Double
shrinkCov sample coef target = scale (1 - coef) sample `add` scale coef target

sampleCov :: Matrix Double -> Herm Double
sampleCov = snd . meanCov

avgVarCov :: Herm Double -> Herm Double
avgVarCov (unSym -> s) = trustSym $ scale s' (ident p)
  where
    p = cols s
    s' = mean $ takeDiag s

naiveCov :: Matrix Double -> Herm Double
naiveCov = avgVarCov . sampleCov

pinCoef :: Double -> Double
pinCoef = max 0 . min 1

lwShrinkageCoef :: Matrix Double -> Double
lwShrinkageCoef x = numer / denom
  where
    n = rows x
    p = cols x
    s = unSym $ sampleCov x
    trs2 = trace $ s <> s
    tr2s = sq $ trace s
    frob ix = let c = x ? [ix] in sq . norm_Frob $ unSym (mTm c) - s
    numer = sumVector kbn (VU.generate n frob) / fromIntegral (sq n)
    denom = trs2 - tr2s / fromIntegral p

rblwShrinkCovCoef :: Int -> Herm Double -> Double
rblwShrinkCovCoef (fromIntegral -> n) (unSym -> s) = numer / denom
  where
    p = cols s
    trs2 = trace $ s <> s
    tr2s = sq $ trace s
    numer = (n - 2) / n * trs2 + tr2s
    denom = (n + 2) * (trs2 - tr2s / fromIntegral p)

rblwShrinkageCoef :: Matrix Double -> Double
rblwShrinkageCoef x = rblwShrinkCovCoef n (sampleCov x)
  where n = rows x

data OASIter = OASIter Double (Herm Double)

data OAS = OAS Int Int (Herm Double) (Herm Double)

oracleStep :: OAS -> Herm Double -> OASIter
oracleStep (OAS n p s f) sj = OASIter pj' sj'
  where
    n' = fromIntegral n
    p' = fromIntegral p
    trss = trace $ unSym sj <> unSym s
    tr2s = sq . trace $ unSym sj
    numer = (1 - 2 / p') * trss  + tr2s
    denom = (n' + 1 - 2 / p') * trss + (1 - n' / p') * tr2s
    pj' = numer / denom
    sj' = scale (1 - pj') s `add` scale pj' f

oracleIter :: Double -> OAS -> Herm Double -> Double
oracleIter tol oas = go 1 . oracleStep oas
  where
    go p (OASIter p' sj')
      | abs (p - p') < tol = p'
      | otherwise = go p' $ oracleStep oas sj'

oracleShrinkCovCoef :: Int -> Herm Double -> Double
oracleShrinkCovCoef n s = oracleIter 1e-6 (OAS n p s f) s
  where
    p = cols $ unSym s
    f = avgVarCov s

oracleShrinkageCoef :: Matrix Double -> Double
oracleShrinkageCoef x = oracleShrinkCovCoef n (sampleCov x)
  where n = rows x

sq :: Num a => a -> a
sq a = a * a

trace :: Matrix Double -> Double
trace = sumVector kbn . takeDiag
