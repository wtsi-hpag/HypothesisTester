# Hypothesis Testing Engine

This is a header-only C++ library which implements a Bayesian odds-ratio computation to comparedifferent statistical models on the same set of data. 

Since statistical models often involve products of very small numbers across high dimensional space (resulting in high susceptibility to floating point errors), we offer some variations on standard integrators to provide more robust estimators of the associated hypervolumes, and hence more accurate estimates of the most appropriate model.