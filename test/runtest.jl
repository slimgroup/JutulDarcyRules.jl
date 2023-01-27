using Test
using Jutul
using JutulDarcyAD
using Flux
using MultiComponentFlash
using Printf
using Random
using LinearAlgebra

include("test_utils.jl")

include("test_conversion.jl")
include("test_jutulState.jl")
include("test_gradient.jl")