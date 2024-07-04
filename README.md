# r-lasla
This repository implements *LASLA*: Locally Adaptive Algorithms for Multiple Testing with Network Structure, with Application to Genome-Wide Association Studies. 
LASLA provides a principled and generic framework that is capable of incorporating network data or multiple samples of auxiliary data from related source domains; 
possibly in different dimensions/structures and from diverse populations. 
LASLA employs a $p$-value weighting approach, utilizing structural insights to assign data-driven weights to individual test points. 
Theoretical analysis shows that LASLA can control FDR with independent or weakly dependent primary statistics, and achieve higher power when the network data is informative.
Efficiency again of LASLA is illustrated through various synthetic experiments and an application to T2D-associated SNP identification.

## Contents

 - `methods/` R codes implementing our methods, benchmarks and some utility functions.
 - `experiments/` Codes to replicate all numerical experiments and corresponding figures discussed in the accompanying paper.
  
  
    
## Prerequisites

Prerequisites libraries:
 - tidyverse
 - progress
