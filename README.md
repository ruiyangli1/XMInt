# XMInt

Model Selection for Multivariate Mediation with Interaction


## About 

A mediation model examines how an independent variable or an exposure (X) affects a dependent variable (Y) through one or more intervening variables or mediators (M). 

This package aims to identify the mediators and the exposure by mediator interactions in the multivariate mediators setting while addressing the underlying hierarchical structure between the main effect and interaction.


## Installation

Users may need to install gfortran before installing our package. Detailed installation instructions can be found at <https://gcc.gnu.org/wiki/GFortranBinaries>. 


## Example

```
## install package
devtools::install_github("ruiyangli1/XMInt")

## load package
library(XMInt)

## data generation
data = dat_gen(N = 200, V = 100, es = 1)
X = data$X; Y = data$Y; M = data$M

## select model
result = XMInt_select(X,Y,M)

## results: selected mediator(s)
result$selected_mediator
## results: selected interaction(s)
result$selected_interaction 
```
