# XMInt

Model Selection for Exposure-Mediator Interaction


## About 

A mediation model examines how an independent variable or an exposure (X) affects a dependent variable (Y) through one or more intervening variables or mediators (M). 

This package conducts the mediation analysis to identify the mediators and the exposure by mediator interactions in the high-dimensional mediators setting while addressing the underlying hierarchical structure between the main effects and interactions.


## Installation

Users may need to install gfortran before installing our package. Detailed installation instructions for gfortran can be found at <https://gcc.gnu.org/wiki/GFortranBinaries>. 

```
## install package
# install.packages("devtools")
devtools::install_github("ruiyangli1/XMInt")

## load package
library(XMInt)
```


## Usage

```
## data generation: true mediators are M1,M2,M3, true interaction is X*M1
data = dat_gen(N = 200, V = 50, es = 1, seed = 1)
X = data$X; Y = data$Y; M = data$M

## select model
result = XMInt_select(X,Y,M)

## selected mediator(s)
result$selected_mediator

## selected interaction(s)
result$selected_interaction 
```
