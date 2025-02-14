#!/usr/bin/bash

CLEAN=(0 1)
PRIORS=("Normal(100.0, 25.0)" "Uniform(50, 300)" "RefOffExponential(50.0, 50.0, 1)")

parallel julia origin_barracudas.jl ::: ${CLEAN[@]} ::: "${PRIORS[@]}"
