OptimizationServices.m
======================

Prototype nonlinear Matlab interface to generate OSiL-format optimization problem instances. Currently only handles scalar variables and can only determine linear coefficients for simple constraints, but appears to work at least for small problems. Based on simple operator overloading and fairly thin wrappers around the OSiL XML content. Problem definition method is similar to Yalmip, for every new variable you get an object that you can perform standard mathematical operations with.
