# TransformUtils.jl


[![TransformUtils](http://pkg.julialang.org/badges/TransformUtils_0.4.svg)](http://pkg.julialang.org/?pkg=TransformUtils&ver=0.4)  [![codecov.io](https://codecov.io/github/dehann/TransformUtils.jl/coverage.svg?branch=master)](https://codecov.io/github/dehann/TransformUtils.jl?branch=master)

[![TransformUtils](http://pkg.julialang.org/badges/TransformUtils_0.5.svg)](http://pkg.julialang.org/?pkg=TransformUtils&ver=0.5)


Lie groups and algebra with some quaternions

## Introduction

This package is a growing collection of Lie Group/Algebra, Quaternion, Euler, AxisAngle representations. Including convert functions between each, and overloading operators for those sets. The package already includes exponential and logarithm maps for Lie SO(3). SE(2) and SE(3) being written as we go.

### Interesting usage

Supports mangle products, for example (using identity constructors):

    julia> pp = SO3(0) * Quaternion(0) * so3(0.1*randn(3)) * AngleAxis(0)

Or maybe you want to compare that we still have identity rotation

    julia> compare(SO3(0), pp) # returns true or false

## Install

    Pkg.add("TransformUtils")

# To do's

Mangled compare functions

    compare(Quaternion(0), SO3(0))
