# Bilevel Toolbox
## Introduction
This is a MATLAB Toolbox designed to test different bilevel optimization problems of the form

---
## Lower Level Problem
In order to define the lower level problem, we need to create a struct that contains the following methods: SOLVE and EVAL

---
## Upper Level Problem
This struct contains the ADJOINT method as well as a EVAL that uses th eoutput from the lower level problem.

---
## Bilevel Solver
Once both the upper and lower level problems have been properly defined, we can run the bilevel solver. To call this solver some previous parameter configurations are needed.
