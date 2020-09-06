Robust Irregular Tensor Factorization and Completion for Temporal Health Data Analysis

By Yifei Ren, Jian Lou, Li Xiong, and Joyce C. Ho. 2020. Code for ``Robust Irregular TensorFactorization and Completion for Temporal Health Data Analysis". In The 29th ACM International Conference on Information and Knowledge Management (CIKM ’20), 2020.

This repository is designed for a robust PARAFAC2 tensor factorization method
for irregular tensors with a new low-rank regularization function to handle potentially missing and erroneous entries in the input tensor.

Before running the codes you need to import Tensor Toolbox Version 2.6 which can be downloaded from: https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html

To start, you need to run: "Run_This.m" file. 
You can select if you want to use the robust version or not:
Smooth_COPA: Smooth PARAFAC2 where smoothness apply to U_k factor matrix.
Robust_Smooth_COPA : Our robust Repair model.


Here is the lists of Repair functions:

 * calculate_fit: -Compute the fit for PARAFAC2 tensor input.

 * claculate_norm: -Compute the norm of a PARAFAC2 tensor.
 
  * claculate_norm_observe:  -Compute the norm of a PARAFAC2 tensor with error & missing entries.

 * MSplineBasis:  -This function produce the spline function for subject X_k.

 * Robust_Smooth_COPA: Robust Smooth PARAFAC2 where smoothness apply to U_k factor matrix.

 * Robust_fastADMM:  -Compute the admm for each mode of a tensor.

 * Robust_COPA_optimizer  -This function is designed to  optimize H,W (S_k) and V.

If you find any bug or error in the codes please send an email to: yifei.ren2@emory.edu
