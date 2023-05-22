# Create a class that stores a stochastic volatility model including errors

setClass("strucvolmodel", representation(pars = "numeric", errors = "numeric",
                                         details = "list"))