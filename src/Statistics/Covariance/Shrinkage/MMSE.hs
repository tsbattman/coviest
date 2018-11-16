
module Statistics.Covariance.Shrinkage.MMSE (
    mmseCov
  , shrinkCov
  , sampleCov
  , naiveCov
  , lwShrinkageCoef
  , rblwShrinkageCoef
  , oracleShrinkageCoef
  ) where

import Numeric.LinearAlgebra hiding ((<>))
import Numeric.Sum (kbn, sumVector)
import Statistics.Sample (mean)
import qualified Data.Vector.Unboxed as VU

mmseCov :: (Matrix Double -> Double) -> Matrix Double -> Herm Double
mmseCov shrink x = shrinkCov s (min 1 p) f
  where
    s = sampleCov x
    f = naiveCov x
    p = shrink x

shrinkCov :: Herm Double -> Double -> Herm Double -> Herm Double
shrinkCov sample coef target = scale (1 - coef) sample `add` scale coef target

sampleCov :: Matrix Double -> Herm Double
sampleCov = snd . meanCov

naiveCov :: Matrix Double -> Herm Double
naiveCov x = trustSym . scale s' $ ident p
  where
    p = cols x
    s' = mean . takeDiag . unSym $ sampleCov x

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

rblwShrinkageCoef :: Matrix Double -> Double
rblwShrinkageCoef x = numer / denom
  where
    n = fromIntegral $ rows x
    p = cols x
    s = unSym $ sampleCov x
    trs2 = trace $ s <> s
    tr2s = sq $ trace s
    numer = (n - 2) / n * trs2 + tr2s
    denom = (n + 2) * (trs2 - tr2s / fromIntegral p)

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

oracleShrinkageCoef :: Matrix Double -> Double
oracleShrinkageCoef x = oracleIter 1e-6 (OAS n p s f) s
  where
    n = rows x
    p = cols x
    s = sampleCov x
    f = trustSym $ scale (mean (takeDiag (unSym s))) (ident p)

sq :: Num a => a -> a
sq a = a * a

trace :: Matrix Double -> Double
trace = sumVector kbn . takeDiag
