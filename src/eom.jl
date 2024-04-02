module GoldbergerWiseEoM
module Consts 
    export yₘ, u, ϕP, ϕT, k, M_IR, γ²₀
    # parameters
    ### [m]=-1
    const yₘ = 1e0π #overall normalization s.t. y_m * M_Pl = pi
    
    ### [m]=0
    # l² = kappa^2 * phiP^2 / 2 reflects the strength of backreaction
    
    ### [m]=1
    const M_Pl = 1.0 #Plank mass
    const u  = 1.0e-1 #  u = log(ϕT / ϕP)/yₘ, this parameter should be fine-tuned to satisfy k*ym ~ O(50), but this causes instability by inrtoducing such a big hierarchy in numerical computation. Therefore here k*ym is set at O(10)
    const k = 37u #pp13 below eq(6.6)
    const M_IR = exp(-k*yₘ) * M_Pl #IR brane scale;(with M_Pl=1)
    const γ²₀ = 4k+2u
    # γ² initially is at large gamma limit
    
    ### [m]=3/2
    const ϕP = 1.e-1 # The scalar field value at Plank-brane
    const ϕT = exp(-u * yₘ)*ϕP
    

end

module AffliatedFunctions
    using ..Consts: yₘ, u, ϕP, ϕT
    export ϕ0, A, A′, A′′, V, V′, V′′, λP′, λT′, λP′′, λT′′, W, W′, α
    #static profile
    @inline ϕ0(  y) = ϕP * (ϕT/ϕP)^(y/yₘ) #ϕ0' = -u ϕ0
    @inline A(   y, l², k, γ²) = k * y + l²/6 * (ϕT/ϕP)^(2y/yₘ)
    @inline A′(  y, l², k, γ²) = k     + l²/6 * (ϕT/ϕP)^(2y/yₘ)* (-2u)
    @inline A′′( y, l², k, γ²) =         l²/6 * (ϕT/ϕP)^(2y/yₘ)* 4u^2
    @inline V′(  y, l², k, γ²) = u*(4k + u)*ϕ0(y) - 2/3*u^2 * 2l²/ϕP^2*ϕ0(y)^3 #κ^2 = 2l²/ϕP^2
    @inline V′′( y, l², k, γ²) = u*(4k + u -2u*2l²*(ϕT/ϕP)^(2y/yₘ))
    @inline λP′( φ, l², k, γ²) = -2u * (ϕP + φ)
    # @inline λP′( φ, l², k, γ²) = -2u * ϕP  + 2γ² * φ # no varphi needed, ϕ≠ϕₚ
    # @inline λP′( φ, l², k, γ²) = -2u * ϕP * φ + 2γ² * φ
    @inline λT′( φ, l², k, γ²) =  2u * (ϕT + φ)
    # @inline λT′( φ, l², k, γ²) =  2u * ϕT  + 2γ² * φ # no varphi needed, ϕ≠ϕₜ
    # @inline λT′( φ, l², k, γ²) =  2u * ϕT * φ + 2γ² * φ
    @inline λP′′(φ, l², k, γ²) = 2γ²
    @inline λT′′(φ, l², k, γ²) = 2γ²

    # superpotential
    W(ϕ, κ, k, u) = 6k/κ^2  - u*ϕ^2
    W′(ϕ, κ, k, u) =        - 2u*ϕ
    V(W, ϕ, κ, k, u) = 1/8 * (W′(ϕ, κ, k, u))^2 - κ^2 / 6 * W(ϕ, κ, k, u)^2

    # for convenience
    α(l²) = 3/2/(l²*u)
end # module AffliatedVariables

module EoMs
    using ..Consts: u, ϕP, yₘ
    using ..AffliatedFunctions: A, A′, A′′, α, ϕ0
    export P, Q, φ, φ′#, RR, SS, getφ′T, getφ′P
    # F'' + P F' + Q F = 0
    """
        P
    P for the EoM: F'' + P F' + Q F = 0
    """
    P(y, l², k, γ², m²) = 2u - 2A′(y, l², k, γ²)
    P(y, params) = P(y, params...)
    """
        Q
    Q for the EoM: F'' + P F' + Q F = 0
    """
    Q(y, l², k, γ², m²) = m² * exp(2A(y, l², k, γ²)) - 4A′′(y, l², k, γ²) - 4u*A′(y, l², k, γ²)
    Q(y, params) = Q(y, params...)
   

    using ForwardDiff
    function φ(F′, F, y, l², k, γ², m²)
        return -3ϕP^2/(2u*l²*ϕ0(y)) * (F′ - 2A′(y, l², k, γ²)*F)  #(3.12)
    end
    φ(F′, F, l², k, γ², m²) = y->φ(F′, F, y, l², k, γ², m²)
    φ(F′, F, params) = φ(F′, F, params...)
    φ′(F′, F, l², k, γ², m²) = y->ForwardDiff.derivative(φ(F′, F, l², k, γ², m²), y)
    φ′(F′, F, params) = φ′(F′, F, params...)
end # module EoMs

module BCs
    using ..Consts: u, ϕP, ϕT, yₘ
    using ..AffliatedFunctions
    using ..EoMs: φ, φ′ #getφ′T, getφ′P, 
    export dφT, Δφ′Ttest, Δφ′Ptest
    #BCs
    # non-perturbative gamma:
    """
        dφP
    given the value of F and F′ at the Plank brane, return φ′ wrt to the BC at Plank brane
    """
    function dφP(F′, F, l², k, γ², m²)
        φ = φP(F′, F, l², k, γ²)
        return 0.5λP′′(φ, l², k, γ²) * φ - 2u * ϕP * F  #(3.14)
    end
    dφP(F′, F, params) = dφP(F′, F, params...)
    φP(F′P, FP, l², k, γ²)  = -3ϕP^2/(2u*l²*ϕ0(0)) * (F′P - 2A′(0, l², k, γ²)*FP)  #(3.12)
    """
        Δφ′P
    given the value of F and F′ at the Plank brane, return the error of φ′to satisfy the BC at Plank brane
    """

    """
    dφT
    given the value of F and F′ at the TeV brane, return φ′ satisfying the BC at TeV brane
    """
    function dφT(F′, F, l², k, γ², m²) 
        φ = φT(F′, F, l², k, γ²)
        return -0.5λT′′(φ, l², k, γ²) * φ - 2u * ϕT * F  #(3.14)
    end
    dφT(F′, F, params) = dφT(F′, F, params...)
    φT(F′T, FT, l², k, γ²)  = -3ϕP^2/(2u*l²*ϕ0(yₘ)) * (F′T - 2A′(yₘ, l², k, γ²)*FT)  #(3.12)

    """
        Δφ′T
    given the value of F and F′ at the TeV brane, return the error of φ′to satisfy the BC at TeV brane
    """
    function Δφ′Ttest(F′, F, params)
        φ′T = dφT(F′, F, params...)
        return φ′T - φ′(F′, F, params)(yₘ)        
    end
    function Δφ′Ptest(F′, F, params)
        φ′P = dφP(F′, F, params...)
        return φ′P - φ′(F′, F, params)(0)        
    end

    """
        dFP
    given the value of F at the Plank brane, return F′ satisfying the BC at Plank brane
    """
    # dFP(FP, l², k, γ², m²) = a(l², k, γ², m²)*FP
    # dFP(FP, params) = dFP(FP, params...)


end # module BCs

end # module GoldbergerWiseEoM
