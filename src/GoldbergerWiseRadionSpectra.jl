module GoldbergerWiseRadionSpectra
export extract_sol, search, m2s_analytical

using DifferentialEquations: ODESolution, ODEProblem, SecondOrderODEProblem, solve, ImplicitEuler
# using Plots
# using LaTeXStrings
include("eom.jl")
using .GoldbergerWiseEoM.EoMs
import .GoldbergerWiseEoM.Consts: yₘ, k, M_IR, γ²₀, u
using .GoldbergerWiseEoM.AffliatedFunctions: ϕ0, A, A′
using .GoldbergerWiseEoM.BCs

function eom!(ddf, df, f, params, y)
    F  = f[1]
    F′ = df[1]
    dF′ = -P(y,params)*F′ -Q(y,params)*F#(3.17)
    ddf[1]=dF′
end

"""
    solveODE_LR(F0, params, y0)
# Arguments
- `F0`: initial value of F′ and f
- `params`: parameters of the model(l², k, γ², m²)
- `y0`: the initial value of y
"""
function solveODE_LR(F0, params, y0)
    F′, F = F0
    yspanL = (y0, 0)
    yspanR = (y0, yₘ)

    probL = SecondOrderODEProblem(eom!,[F′], [F],yspanL, params)
    probR = SecondOrderODEProblem(eom!,[F′], [F],yspanR, params)
    
    solL = solve(probL, ImplicitEuler())
    solR = solve(probR, ImplicitEuler())
    sol(y) = y < y0 ? solL(y) : solR(y)

    return sol
end
function extract_sol(props)
    l2, g2, m2, dm2, F′ = props
    params = (l2, k, g2*γ²₀, m2*M_IR^2)
    yspan = (yₘ, 0)
    # return solve(SecondOrderODEProblem(eom!,[F′], [1.],yspan, params), ImplicitEuler())
    return solveODE_LR([F′, 1.], params, 0.99yₘ)
end
function errBCwithφ(params; F0′, F0=1., y0 )
    Fsol = solveODE_LR([F0′, F0], params, y0)
    FT′,FT = Fsol(yₘ)
    FP′,FP = Fsol(0)
    errR = Δφ′Ttest(FT′, FT, params)
    errL = Δφ′Ptest(FP′, FP, params)
    return errL, errR
end
function searchYM(l2, g2, y0, k=k)
    (Fp, m2) -> errBCwithφ((l2, k, g2*γ²₀, m2*M_IR^2); F0′=Fp, F0=1., y0=y0)
end
function m2s_analytical(l2, g2)
    return m2s = 4l2*(2k+u)*u^2/(3k)*(1-exp(2k*yₘ))/(1-exp((4k+2u)*yₘ)) * (1 .- 1 /g2) / M_IR^2
end
include("search_alg.jl")
function search(l2=1e-3, g2=1e2; y0=0.99yₘ, niters=5, scale0 = 10, m2scale = 1e-4, Fprange = (2., 13.), tol=1e-6, nresol=10)
    # both Fp and m2 leaner scaled
    prob = searchYM(l2, g2, y0, k)
    n = nresol
    m20 = 0
    Fp1 = bisection(Fp->prob(Fp, m20)|>first, Fprange...)
    Fp2 = bisection(Fp->prob(Fp, m20)|>last, Fprange...)
    Fp0 = 0.5(Fp1+Fp2)
    Fpscale = 0.8abs(Fp1-Fp2)
    # extract the hypersurface constrained by BCs
    for i in 1:niters
        Fp, m2, err = sample(prob, Fp0, Fpscale, m20, m2scale, n)
        Fp0_old, m20_old = Fp0, m20
        Fp0, m20, hypersurfaces = fit(Fp, m2, err)
        if isnan(m20) || isnan(Fp0)  
            scale=1/scale0
            Fp0, m20 = Fp0_old, m20_old 
        else
            scale=scale0
        end
        m2scale, Fpscale = rescale(m2scale, Fpscale, m20, Fp0, hypersurfaces, scale)
        m2scale, Fpscale = max(m2scale, abs(m20-m20_old)), max(Fpscale, abs(Fp0-Fp0_old)) # avoid to jump too far from the previous region
        tol > m2scale && break
    end
    Fp, m2, err = sample(prob, Fp0, Fpscale, m20, m2scale, n)
    Fp0, m20, hypersurfaces = fit(Fp, m2, err)
    return m20, m2scale/n, Fp0
end

end
