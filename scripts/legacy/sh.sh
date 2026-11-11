#!/bin/sh 
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random/intersection_cuts_relaxed_y/complete.batch > log_campinas_random_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random-more-components/intersection_cuts_relaxed_y/complete.batch > log_campinas_random_more_components_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random-unitary-profits/intersection_cuts_relaxed_y/complete.batch > log_campinas_random_unitary_profits_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random-more-components-unitary-profits/intersection_cuts_relaxed_y/complete.batch > log_campinas_random_more_components_unitary_profits_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse/intersection_cuts_relaxed_y/complete.batch > log_campinas_sparse_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse-more-components/intersection_cuts_relaxed_y/complete.batch > log_campinas_sparse_more_components_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse-unitary-profits/intersection_cuts_relaxed_y/complete.batch > log_campinas_sparse_unitary_profits_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse-more-components-unitary-profits/intersection_cuts_relaxed_y/complete.batch > log_campinas_sparse_more_components_unitary_profits_complete
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random/brkga.batch > log_campinas_random_brkga
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random-more-components/brkga.batch > log_campinas_random_more_components_complete_brkga
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random-unitary-profits/brkga.batch > log_campinas_random_unitary_profits_complete_brkga
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-random-more-components-unitary-profits/brkga.batch > log_campinas_random_more_components_unitary_profits_brkga
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse/brkga.batch > log_campinas_sparse_brkga
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse-more-components/brkga.batch > log_campinas_sparse_more_components_brkga
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse-unitary-profits/brkga.batch > log_campinas_sparse_unitary_profits_brkga
julia --threads=auto --project=. src/run.jl --batch batchs/campinas-sparse-more-components-unitary-profits/brkga.batch > log_campinas_sparse_more_components_unitary_profits_brkga
