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

"""
    extract_sol takes the properties of the radion profile and returns the solution of the radion profile.
# Arguments
- `props`: the properties of the radion profile, including l², g², m², dm², F′
# Returns
- `sol`: the solution of the radion profile
# Example
```julia
l2 = 1e-3
g2 = 1e3
extract_sol([l2, g2, search(l2, g2)...])
```
"""
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

"""
    m2s_analytical calculates the analytical (perturbative) value of the radion mass squared.
# Arguments
- `l2`: the backreaction parameter
- `g2`: the quadratic brane potential parameter, in unit of the critical value g^2_0 = 4k+2u
# Returns
- `m2s`: the analytical value of the radion mass squared, in unit of M_IR
# Example
```julia
m2s_analytical(1e-3, 1e3)
```
"""
function m2s_analytical(l2, g2)
    return m2s = 4l2*(2k+u)*u^2/(3k)*(1-exp(2k*yₘ))/(1-exp((4k+2u)*yₘ)) * (1 .- 1 /g2) / M_IR^2
end
include("search_alg.jl")

"""
    find the eigenmass of the radion (its error bound and the corresponding Fp) by searching the parameter space.
# Arguments
- `l2`: the backreaction parameter
- `g2`: the quadratic brane potential parameter, in unit of the critical value g^2_0 = 4k+2u
- `y0`: the initial shooting value of y position
- `niters`: the maximum number of iterations
- `scale0`: the default scale to zoom in/out the patch of the parameter space
- `m2scale`: the initial patch size of m²
- `Fprange`: the default range searching for Fp that satisfies boundary conditions
- `tol`: the tolerance of the patch size of m²
- `nresol`: the number of points to sample each dimension of the patch of parameter space
# Returns
- `m20`: the eigenmass of the radion, in unit of M_IR
- `m2scale/n`: the error bound of the eigenmass
- `Fp0`: the corresponding Fp
# Example
```julia
m20, m2err, Fp0 = search(1e-3, 1e3)
```
# Note
When the bulk equation is solved without using any of boundary conditions, there are several integration constants and the square mass eigenvalue left to be determined.
In this sense, we say that the solutions of the bulk equation are parameterized by squared mass eigenvalue m^2 and these integration constants C_0, D_0^1, etc..
In other words, each point in this original parameter space of (m^2, C_0, D_0^1, ...) can represent a solution of the bulk equation.

But the "true" solution we really want was a small subset inside the original parameter space, which in turn should be constrained by the boundary conditions.

Nevertheless, the boundary conditions can be very complicated to deal with, e.g. they locate in different locations of the bulk(i.e. become boundary value problem rather than the initial value problem), they mix the field value and its derivative and even other fields.

So how can we translate the boundary conditions into constraints in the parameter space explicitly?
Since we do not have an analytical form of the solution yet, we cannot directly specify a point in the original parameter space. So we need a more suitable description(in the sense of coordinates transformation) of this parameter space, which can be easily and directly controlled in numerical simulations.
Here we choose to introduce an "intermediate" set of boundary conditions BC0( F(y_0), F'(y_0)), which directly specifies the field value F(y_0) and its derivative F'(y_0) at y=y_0, where y_0 now is an arbitrary position chosen in [0,\pi R]. By this step, we can obtain a useful effective description--a practical parameter space, by squared mass m^2 and the first order derivative of field F'(y_0). Here are the reasons why it's valid and useful:

Firstly, we can verify that only F(y_0) and F'(y_0) should be enough to determine all integration constants(when squared mass is also specified), since the bulk equation is a 2nd order ODE. So each element in the "intermediate" boundary condition set with a specified squared mass value maps to a point in the original parameter space  (m^2, C_0, D_0^1, ...). Therefore, we have an effective description of (m^2, C_0, D_0^1, ...), by the practical one (m^2, F(y_0), F'(y_0)). As a sidenote, F(y_0) is actually redundant because it corresponds to the overall normalization, so we conclude that only (m^2, F'(y_0)) should be enough to parameterize all solutions of the bulk equation.

Secondly, there is a huge advantage of taking such a "intermediate" set of boundary conditions, because it is in the simplist form for numerical simulation--the initial value problem, which can be solved with many powerful numerical solvers for very complicated differential equations.

Thirdly, the freedom of setting y_0 help us to balance the numerical error when the solution have exponential behavior. For example, if the solution F( y) include exp(ky) term where k is a positive constant, and we try to have a initial guess at y_0 = 0, then the numerical error from F(y_0) and F'(y_0) would be exponentially amplified when evolving through (0,\pi R). On the other hand, if we set y_0=\pi R, then the numerical error would be suppressed, when evolving from \pi R to 0.

Now, we have constructed the practical parameter space (m^2, F'(y_0)) that can represent all possible solutions of the bulk equation. All parameters m^2 and F'(y_0) can be easily manipulated in numerical simulations.
We continue to answer the question how can we translate the boundary conditions into constraints in the parameter space explicitly.
The strategy in numerical analysis here is nothing but "TRY AND ERROR", in a little bit tricky and efficient way.
Each time we choose a specific point in the practical parameter space of (m^2, F'(y_0)), we can rebuild the corresponding solution of the bulk equation, and see at boundaries whether its field value F and derivative F' satisfies the boundary condition BC(F,F')=0.
For each boundary condition, the points satisfying it form a curve in the parameter space. And the "true" solution should locate in the intersection of all curves, by means that the "true" solution should satisfy all boundary conditions simutaneously.
So to get the "true" solution, we only need to find all curves in the practical parameter space (m^2, F'(y_0)) corresponding to the boundary conditions at Planck and the TeV branes, then extract their intersections.
But more practically, to reduce the cost of computation, we do not scan the whole practical parameter space of (m^2, F'(y_0)).
Instead, we scatter some initial guesses on a small patch (but large enough to cover at least one piece of each curve corresponding to boundary conditions) of the practical parameter space. We would like to zoom in/out the patch to let it focus on the intersection point that represents the "true" solution. This can be realized in the following way:

1). Fit for each boundary condition in the practical parameter space and find the intersection of fitted curves.

2). If the intersection point is outside the current patch, zoom out to cover it. Otherwise, zoom in to focus on it.

3). Repeat step 1 and 2 until the patch is zoomed in so that its size is smaller than the numerical error tolerance.

The output of this algorithm is a small enough patch that covers the intersection point of curves corresponding to boundary conditions, which represents the "true" solution satisfying all boundary conditions simutaneously. By setting the numerical error tolerance, we can locate the intersection point of all curves of boundary conditions at any precision.
"""
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
