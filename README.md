# ClinicalTrialOptm.jl

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://Hua-Zhou.github.io/ClinicalTrialOptm.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://Hua-Zhou.github.io/ClinicalTrialOptm.jl/dev) | [![build Actions Status](https://Hua-Zhou.github.io/ClinicalTrialOptm.jl/workflows/CI/badge.svg)](https://Hua-Zhou.github.io/ClinicalTrialOptm.jl/actions)  | [![Coverage Status](https://coveralls.io/repos/github/Hua-Zhou/ClinicalTrialOptm.jl/badge.svg?branch=main)](https://coveralls.io/github/Hua-Zhou/ClinicalTrialOptm.jl?branch=main) [![codecov](https://codecov.io/gh/Hua-Zhou/ClinicalTrialOptm.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Hua-Zhou/ClinicalTrialOptm.jl) |  


## The Problem

<img alt="clinicaltrialimage" src="https://github.com/Hua-Zhou/ClinicalTrialOptm.jl/assets/110351617/bfebe101-df73-43fc-af0b-a3d41525586f">

In the last five years, clinical trials have become increasingly more difficult to conduct due to staffing, budget, and protocol complications. According to the 2020 Tufts University Impact Report, more than 20% of clinical trials fail to recruit enough patients in time. Biomedical and phameceutical companies rely on trial data to progress in treatment development, making effective clinical trial design necessary to address.

ClinicalTrialOptm.jl solves these multi-center, multi-state recruitment problems using mixed-integer algorithms. It seeks to optimize the number of clinics for each country, minimizing cost while maintaining high probabilities of successful recruitment. Details on the calculations are described in the paper: (Insert paper here)

## Installation 
ClinicalTrialOptm.jl supports Julia v1.7 or later. See documentation for usage. It is not yet registered and can be installed, in the Julia Pkg mode, by the following command:

```{julia}
(@v1.8) Pkg> add https://github.com/Hua-Zhou/ClinicalTrialOptm.jl.git
```
