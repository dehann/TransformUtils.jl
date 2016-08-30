# TransformUtils.jl


[![TransformUtils](http://pkg.julialang.org/badges/TransformUtils_0.4.svg)](http://pkg.julialang.org/?pkg=TransformUtils&ver=0.4)  [![codecov.io](https://codecov.io/github/dehann/TransformUtils.jl/coverage.svg?branch=master)](https://codecov.io/github/dehann/TransformUtils.jl?branch=master)

[![TransformUtils](http://pkg.julialang.org/badges/TransformUtils_0.5.svg)](http://pkg.julialang.org/?pkg=TransformUtils&ver=0.5)


Lie groups and algebra with some quaternions

## Introduction

This package is a growing collection of Lie Group/Algebra, Quaternion, Euler, AxisAngle representations. Including convert functions between each, and overloading operators for those sets. The package already includes exponential and logarithm maps for Lie SO(3). SE(2) and SE(3) being written as we go.

Tests are still pending, but the code has been used quite a lot in previous projects etc.

## Install

    Pkg.add("TransformUtils")
