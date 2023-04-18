using Test
using Printf
using LinearAlgebra
mean(x) = sum(x)/length(x)

function log_division(a, b)
    if b == 0
        return a == 0 ? NaN : Inf
    end
    return log(a / b)
end

function grad_test(f, x0, dx, dJdx; dJ=nothing, maxiter=6, h0=5e-2, stol=1e-1, hfactor=8e-1)
    if !xor(isnothing(dJdx), isnothing(dJ))
        error("Must specify either dJdx or dJ") 
    end
    if isnothing(dJ)
        dJ = dot(dJdx, dx)
    end
    J0 = f(x0)
    h = h0

    log_factor = log(hfactor)
    expected_f1 = 1e0 / hfactor
    expected_f2 = 1e0 / hfactor ^2e0

    err1 = zeros(Float64, maxiter)
    err2 = zeros(Float64, maxiter)
    Js = zeros(Float64, maxiter)
    @printf("%11s, %12s | %11s, %11s | %11s, %11s | %12s, %12s | %12s %12s %12s \n", "h", "dJ", "e1", "e2", "factor1", "factor2", "rate1", "rate2", "Finite-diff", "T1 approx", "T2 approx")
    @printf("%11s, % 12.5e | %11s, %11s | %11.5e, %11.5e | % 12.5e, % 12.5e | \n", "", dJ, "0", "0", expected_f1, expected_f2, 1, 2)
    @printf("%11s, %12s | %11s, %11s | %11s, %11s | % 12s, % 12s \n", "___________", "___________", "___________", "___________", "___________", "___________", "____________", "____________")
    all_info = []
    for j=1:maxiter
        J = f(x0 + h*dx)
        Js[j] = J
        err1[j] = norm(J - J0, 1)
        err2[j] = norm(J - J0 - h*dJ, 1)
        j == 1 ? prev = 1 : prev = j - 1

        dJ_est = (J - J0) / h
        α = expected_f2
        dJ_est1 = (α*Js[prev] - J + (1-α)*J0) / (h * (α/hfactor - 1))
        α = -expected_f2
        dJ_est2 = (α*Js[prev] - J + (1-α)*J0) / (h * (α/hfactor - 1))

        rate1 = log_division(err1[j], err1[prev]) / log_factor
        rate2 = log_division(err2[j], err2[prev]) / log_factor

        info = (h, h*norm(dJ, 1), err1[j], err2[j], err1[prev]/err1[j], err2[prev]/err2[j], rate1, rate2, dJ_est, dJ_est1, dJ_est2)
        push!(all_info, info)
        @printf("%11.5e, % 12.5e | %11.5e, %11.5e | %11.5e, %11.5e | % 12.5e, % 12.5e | % 12.5e, % 12.5e, % 12.5e \n", info...)
        h = h * hfactor
    end

    println()
    @printf("%11s, %12s | %11s, %11s | %11s, %11s | %12s, %12s | %12s %12s %12s \n", "h", "dJ", "e1", "e2", "factor1", "factor2", "rate1", "rate2", "Finite-diff", "T1 approx", "T2 approx")
    @printf("%11s, % 12.5e | %11s, %11s | %11.5e, %11.5e | % 12.5e, % 12.5e | \n", "", dJ, "0", "0", expected_f1, expected_f2, 1, 2)
    @printf("%11s, %12s | %11s, %11s | %11s, %11s | % 12s, % 12s \n", "___________", "___________", "___________", "___________", "___________", "___________", "____________", "____________")
    for j=1:maxiter
        info = all_info[j]
        @printf("%11.5e, % 12.5e | %11.5e, %11.5e | %11.5e, %11.5e | % 12.5e, % 12.5e | % 12.5e, % 12.5e, % 12.5e \n", info...)
        h = h * hfactor
    end

    rate1 = err1[1:end-1]./err1[2:end]
    rate2 = err2[1:end-1]./err2[2:end]
    @test mean(rate1) ≥ expected_f1 - stol
    @test mean(rate2) ≥ expected_f2 - stol
end

