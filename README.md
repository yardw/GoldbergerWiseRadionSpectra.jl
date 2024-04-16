# GoldbergerWiseRadionSpectra.jl 
<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yardw.github.io/GoldbergerWiseRadionSpectra.jl/stable/)  -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yardw.github.io/GoldbergerWiseRadionSpectra.jl/dev/)  -->
[![Build Status](https://github.com/yardw/GoldbergerWiseRadionSpectra.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yardw/GoldbergerWiseRadionSpectra.jl/actions/workflows/CI.yml?query=branch%3Amain) 
<!-- [![Coverage](https://codecov.io/gh/yardw/GoldbergerWiseRadionSpectra.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yardw/GoldbergerWiseRadionSpectra.jl) -->

# Installation
To install the package, you can clone the repository manually or by running the following command in the terminal.
```bash
git clone https://github.com/yardw/GoldbergerWiseRadionSpectra.jl.git
```
Enter the **root** of this repository, run `julia` (if you have not installed Julia yet, you can get it from [here](https://julialang.org/downloads/)).
In the Julia REPL, you can enter the package mode by pressing `]` key. 
You will see the prompt changes to
```bash
julia> ]
(@v1.*) pkg>
```
where `(@v1.*)` represents your base environment.
Then `activate` the current package environment and `instantiate` the package.
```jldoctest
activate .
instantiate
```
The above commands will automatically install all necessary dependencies and now the package is ready to be used.

# Usage
After instatiation, you can use the package in the Julia REPL or in a Julia script by running the following command.
```julia
using GoldbergerWiseRadionSpectra
```
The package provides following functions:
- `search` : search for the radion mass spectrum in the Goldberger-Wise model, with respect to given parameters
    - inputs:
        - `l2`: the backreaction parameter
        - `g2`: the quadratic brane potential parameter, in unit of the critical value g^2_0 = 4k+2u
    - outputs:
        - `m20`: the squared mass of the radion
        - `m2err`: the error bound of the squared mass
        - `Fp`: the first-order derivative of the radion profile at the initial position. It can be used to rebuild the radion profile.
- `m2s_analytical` : calculate the perturbative mass spectrum of the radion in the Goldberger-Wise model, with respect to `l2` and `g2`
- `extract_sol` : extract the radion profile from the output of `search`
You can also modify other parameters for Goldberger-Wise model, which includes `yₘ, u, ϕP, ϕT, k, M_IR, γ²₀`, in the `eom.jl` file.