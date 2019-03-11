Generic l1 l2 Solver
********************
This solver is a generic lower level solver designed to solve a wide variety of image processing problems. The abstract structure has the following form

\\[ f(x) = \sum_{k=1}^N \alpha_k \|A_k x -z_k\|^2 + \sum_{l=1}^M \beta_l \|B_l x - z_l\| \\]

where \\(\alpha_k,\beta_l\\) are regularization parameters, \\(A_k,B_l\\) are linear operators, \\(N\\) is the number of \\(l2\\) terms to be consider and likewise \\(M\\) is the number of \\(l_1\\) terms.

This abstract structure allows us to formulate image processing problems such as

1. Image Denoising using TV, TGV and ICTV
2. Image Deblurring
3. Image Inpainting

## Solving the abstract model
The numerical strategy chosen for solving this problem is to make use of a *saddle point formulation*. This structure will allow us to experiment with primal-dual based first order algorithms. In particular the implementation will be made using Chambolle-Pock algorithm.

## Saddle Point Formulation
