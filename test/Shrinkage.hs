
module Shrinkage (
    strimmerX
  , testShrinkage
  ) where

import Numeric.LinearAlgebra
import Test.Tasty
import Test.Tasty.HUnit

import Statistics.Covariance.Shrinkage.Strimmer

strimmerX :: Matrix Double
strimmerX = matrix 3 [
    -1.1541680,  0.1229077,  0.2964592
  , -0.4574491, -0.0106610, -0.2497785
  ,  1.2711852,  0.2362810,  0.3653930
  ]

withinTolVec :: Double -> Vector Double -> Vector Double -> Bool
withinTolVec tol x x' = abs (x - x') < scalar tol

withinTolMat :: Double -> Matrix Double -> Matrix Double -> Bool
withinTolMat tol x x' = withinTolVec tol (flatten x) (flatten x')

strimmerpsmallSVD :: Assertion
strimmerpsmallSVD = assertBool "psmallSVD" $ withinTolMat tol u' u
  && withinTolVec tol s' s
  && withinTolMat tol v' v
  where
    tol = 1e-6
    (u, s, v) = psmallSVD strimmerX
    u' = matrix 3 [
         0.6256382, -0.7737046, -0.0997895
      ,  0.2669970,  0.3325596, -0.9044980
      , -0.7330003, -0.5392451, -0.4146388
      ]
    s' = vector [1.7850070, 0.5601167, 0.1101480]
    v' = matrix 3 [
        -0.99495870,  0.09886453,  0.01682235
      , -0.05554307, -0.40358205, -0.91325599
      , -0.08349942, -0.90958635,  0.40703871
      ]

strimmernsmallSVD :: Assertion
strimmernsmallSVD = assertBool "nsmallSVD" $ withinTolMat tol u' u
  && withinTolVec tol s' s
  && withinTolMat tol v' v
  where
    tol = 1e-6
    (u, s, v) = nsmallSVD strimmerX
    u' = matrix 3 [
        -0.6256382,  0.7737046, 0.0997895
      , -0.2669970, -0.3325596, 0.9044980
      ,  0.7330003,  0.5392451, 0.4146388
      ]
    s' = vector [1.7850070, 0.5601167, 0.1101480]
    v' = matrix 3 [
        0.99495870, -0.09886453, -0.01682235
      , 0.05554307,  0.40358205,  0.91325599
      , 0.08349942,  0.90958635, -0.40703871
      ]

strimmerposSVD :: Assertion
strimmerposSVD = assertBool "posSVD" $ withinTolMat tol u' u
  && withinTolVec tol s' s
  && withinTolMat tol v' v
  where
    tol = 1e-6
    (u, s, v) = posSVD strimmerX
    u' = matrix 3 [
        -0.6256382,  0.7737046, 0.0997895
      , -0.2669970, -0.3325596, 0.9044980
      ,  0.7330003,  0.5392451, 0.4146388
      ]
    s' = vector [1.7850070, 0.5601167, 0.1101480]
    v' = matrix 3 [
        0.99495870, -0.09886453, -0.01682235
      , 0.05554307,  0.40358205,  0.91325599
      , 0.08349942,  0.90958635, -0.40703871
      ]

strimmerCovShrink :: Assertion
strimmerCovShrink = assertBool "corpcor cor shrink is wrong" $ withinTolMat 1e-6 cov' cov
  where
    cov = unSym (covShrink strimmerX)
    cov' = matrix 3 [
        1.13853590, 0.06889376, 0.05675350
      , 0.06889376, 0.04389409, 0.03081658
      , 0.05675350, 0.03081658, 0.11359393
      ]

strimmerCorShrink :: Assertion
strimmerCorShrink = assertBool "corpcor cor shrink is wrong" $ withinTolMat 1e-6 cor' cor
  where
    cor = unSym (corShrink strimmerX)
    cor' = matrix 3 [
        1.0000000, 0.3081793, 0.1578126
      , 0.3081793, 1.0000000, 0.4364192
      , 0.1578126, 0.4364192, 1.0000000
      ]

strimmerCorCoef :: Assertion
strimmerCorCoef = assertBool msg $ abs (lambda - coef) < 0.0000001
  where
    lambda = 0.5311773
    coef = corCoef strimmerX
    msg = "expected: " ++ show lambda ++ "\n     got: " ++ show coef

strimmerVarCoef :: Assertion
strimmerVarCoef = assertBool "corpcor var coef is wrong" $ abs (lambda - varCoef strimmerX) < 0.0000001
  where lambda = 0.2910548

strimmerVarShrink :: Assertion
strimmerVarShrink = assertBool "corp var shrink" $ abs (v - varShrink strimmerX) < 0.0000001
  where v = vector [1.13853590, 0.04389409, 0.11359393]

testShrinkage :: TestTree
testShrinkage = testGroup "shrinkage" [
    testCase "strimmer varCoef" strimmerVarCoef
  , testCase "strimmer varShrink" strimmerVarShrink
  , testCase "strimmer corCoef" strimmerCorCoef
  , testCase "strimmer corShrink" strimmerCorShrink
  , testCase "strimmer covShrink" strimmerCovShrink
  , testCase "strimmer psmallSVD" strimmerpsmallSVD
  , testCase "strimmer nsmallSVD" strimmernsmallSVD
  , testCase "strimmer posSVD" strimmerposSVD
  ]
