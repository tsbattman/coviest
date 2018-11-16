
module Main (main) where

import Test.Tasty

import Shrinkage

main :: IO ()
main = defaultMain $ testGroup "all tests" [
    testShrinkage
  ]

