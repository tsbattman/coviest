
module Statistics.Covariance.Shrinkage.Strimmer (
    psmallSVD
  , nsmallSVD
  , posSVD
  , covShrink
  , corShrink
  , corCoef
  , varShrink
  , varCoef
  ) where

import Numeric.LinearAlgebra hiding ((<>))
import Numeric.Sum (kbn, sumVector)
import Statistics.Quantile (quantile, s)
import qualified Data.Vector.Storable as VS

covShrink :: Matrix Double -> Herm Double
covShrink x = trustSym $ (asColumn sc * unSym c) * asRow sc
  where
    sc = sqrt $ varShrink x
    c = corShrink x

psmallSVD :: Matrix Double -> (Matrix Double, Vector Double, Matrix Double)
psmallSVD x = (u', d', v')
  where
    b = unSym $ mTm x
    (_, d, v) = svd b
    tol = 1e-6
    positive = find (> tol) d
    u' = (x <> v') * asRow (recip d')
    d' = sqrt . vector $ map (atIndex d) positive
    v' = v ¿ positive

nsmallSVD :: Matrix Double -> (Matrix Double, Vector Double, Matrix Double)
nsmallSVD x = (u', d', v')
  where
    b = unSym $ mTm (tr x)
    (u, d, _) = svd b
    tol = 1e-6
    positive = find (> tol) d
    u' = u ¿ positive
    d' = sqrt . vector $ map (atIndex d) positive
    v' = (tr x <> u') * asRow (recip d')

posSVD :: Matrix Double -> (Matrix Double, Vector Double, Matrix Double)
posSVD x = (u', d', v')
  where
    (u, d, v) = svd x
    tol = 1e-6
    positive = find (> tol) d
    u' = u ¿ positive
    d' = vector $ map (atIndex d) positive
    v' = v ¿ positive

fastSVD :: Matrix Double -> (Matrix Double, Vector Double, Matrix Double)
fastSVD x
  | n > edgeRatio * p = psmallSVD x
  | edgeRatio * n < p = nsmallSVD x
  | otherwise = posSVD x
  where
    edgeRatio = 2
    n = rows x
    p = cols x

mpower :: Herm Double -> Double -> Herm Double
mpower x p = trustSym $ (evec * asRow (cmap exp $ scale p (cmap log eval))) <> tr evec
  where (eval, evec) = eigSH x

corShrink :: Matrix Double -> Herm Double
corShrink x
  | lambda == 1 || alpha == 0 = trustSym $ ident p
  | alpha == 1 = unitDiag $ scale (1 - lambda) r0
  | lambda == 0 = unitDiag . trustSym $ svdxsV <> unSym (mpower c2 alpha) <> tr svdxsV
  | otherwise = unitDiag . trustSym . scale (exp (alpha * log lambda)) $ ident p - (svdxsV <> f <> tr svdxsV)
  where
    unitDiag :: Herm Double -> Herm Double
    unitDiag powr = trustSym (diag (1 - takeDiag (unSym powr))) `add` powr
    n = rows x
    p = cols x
    alpha = 1 :: Double
    w2 = recip $ fromIntegral n -- w = 1 / n => sum(w^2) =  1 / n
    h1 = 1 / (1 - w2)
    xs = let (u, v) = meanCov x in (x - asRow u) / asRow (sqrt (takeDiag (unSym v)))
    lambda = corCoef x
    r0 = scale h1 $ mTm (xs / sqrt (fromIntegral n))
    (svdxsU, svdxsS, svdxsV) = fastSVD xs
    m = size svdxsS
    utwu = tr svdxsU <> svdxsU / fromIntegral n
    c0 = utwu * asRow svdxsS * asColumn svdxsS
    c1 = scale ((1 - lambda) * h1) c0
    c2 = sym c1
    f = ident m - unSym (mpower (trustSym $ unSym (scale (recip lambda) c2) + ident m) alpha)

corCoef :: Matrix Double -> Double
corCoef x
  | denominator == 0 = 1
  | otherwise = numerator / denominator * h1w2
  where
    n = rows x
    p = cols x
    xs = let (u, v) = meanCov x in (x - asRow u) / asRow (sqrt (takeDiag (unSym v)))
    w2 = recip $ fromIntegral n
    h1w2 = w2 / (1 - w2)
    sw = sqrt . recip $ fromIntegral n
    xsw = scale sw xs
    (xswsvdU, xswsvdS, xswsvdV) = fastSVD xsw
    colSums = vector . map (sumVector kbn) . toColumns
    -- colCumSums = fromColumns . map (VS.postscanl (+) 0) . toColumns
    rowCumSums = fromRows . map (VS.postscanl (+) 0) . toRows
    sE2R = sumVector kbn (flatten (xsw * ((xswsvdU * asRow (cb xswsvdS)) <> tr xswsvdV)))
      - sumVector kbn (sq (colSums (sq xsw)))
    xs2w = scale sw $ sq xs
    sER2 = 2 * (sumVector kbn . flatten
      $ fliprl (takeColumns (p - 1) xs2w) * (rowCumSums . fliprl $ dropColumns 1 xs2w))
    denominator = sE2R
    numerator = sER2 - sE2R

varShrink :: Matrix Double -> Vector Double
varShrink x = scalar (lambda * target) + scale (1 - lambda) v
  where
    lambda = varCoef x
    v = takeDiag . unSym . snd $ meanCov x
    target = quantile s 1 2 v

varCoef :: Matrix Double -> Double
varCoef x
  | denominator == 0 = 1
  | otherwise = max 0 . min 1 $ numerator / denominator * h1w2
  where
    n = rows x
    w2 = recip (fromIntegral n)
    h1 = 1 / (1 - w2)
    h1w2 = w2 / (1 - w2)
    xc = x - asRow (fst (meanCov x))
    v = scale h1 q1
    target = quantile s 1 2 v
    zz = xc * xc
    q1 = fst $ meanCov zz
    q2 = fst (meanCov (zz * zz) ) - q1 * q1
    numerator = sumVector kbn q2
    denominator = sumVector kbn . sq $ q1 - scalar (target / h1)

sq, cb :: Num a => a -> a
sq a = a * a
cb a = a * sq a
