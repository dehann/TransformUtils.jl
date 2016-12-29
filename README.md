# TransformUtils.jl

[![Build Status](https://travis-ci.org/dehann/TransformUtils.jl.svg?branch=master)](https://travis-ci.org/dehann/TransformUtils.jl)
[![codecov.io](https://codecov.io/github/dehann/TransformUtils.jl/coverage.svg?branch=master)](https://codecov.io/github/dehann/TransformUtils.jl?branch=master)
[![TransformUtils](http://pkg.julialang.org/badges/TransformUtils_0.4.svg)](http://pkg.julialang.org/?pkg=TransformUtils&ver=0.4)
[![TransformUtils](http://pkg.julialang.org/badges/TransformUtils_0.5.svg)](http://pkg.julialang.org/?pkg=TransformUtils&ver=0.5)
[![TransformUtils](http://pkg.julialang.org/badges/TransformUtils_0.6.svg)](http://pkg.julialang.org/?pkg=TransformUtils&ver=0.6)

Lie groups and algebra, quaternions, Angle Axis and Euler angles; products and compare also available.

## Introduction

This package is a growing collection of Lie Group/Algebra, Quaternion, Euler, AxisAngle representations. Including convert functions between each, and overloading operators for those sets. The package already includes exponential, logarithm maps and Jacobians for Lie SO(3). SE(3) mostly complete.

### Interesting usage

Supports mangle products, for example (using identity constructors):

    julia> pp = SO3(0) * Quaternion(0) * so3(0.1*randn(3)) * AngleAxis(0)

Or maybe you want to compare against another rotation

    julia> compare(SO3(0), pp) # returns true or false

## Install

    Pkg.add("TransformUtils")

# To do's

Mangled compare functions

    compare(Quaternion(0), SO3(0))

Rework SE(2) to type, and not functions
