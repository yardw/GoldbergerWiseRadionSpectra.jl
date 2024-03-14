# GoldbergerWiseRadionSpectra.jl 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yardw.github.io/GoldbergerWiseRadionSpectra.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yardw.github.io/GoldbergerWiseRadionSpectra.jl/dev/) [![Build Status](https://github.com/yardw/GoldbergerWiseRadionSpectra.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yardw/GoldbergerWiseRadionSpectra.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/yardw/GoldbergerWiseRadionSpectra.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yardw/GoldbergerWiseRadionSpectra.jl)

# Usage
## Installation
To install the package, you can clone the repository manually or by running the following command in the terminal.
```bash
git clone https://github.com/yardw/GoldbergerWiseRadionSpectra.jl.git
```
Enter the **root** of this repository, run `julia` (if you have not installed `julia` yet, you can get it from [here](https://julialang.org/downloads/)).
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