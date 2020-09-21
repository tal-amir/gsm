Code to calculate the Generalized Soft-Min function
===================================================
 
Version 5.1, June 2020
 
Based on the paper:  
Tal Amir, Ronen Basri, Boaz Nadler (2020) - The Trimmed Lasso: Sparse Recovery Guarantees and Practical Optimization by the Generalized Soft-Min Penalty
 
Requirements
------------
**Matlab**  
This program requires Matlab 2018b and onward, but may work with earlier versions.  
Mex binaries are available for Windows and Linux.  
On Mac, the code may be compiled on the target computer, or run avoiding Mex binaries; see below.
 
Usage
-----
Basic usage:  
`>> [mu, theta] = gsm_v5_1(z, k, gamma);`
 
To avoid calling the Mex binary and run purely on Matlab:  
`>> [mu, theta] = gsm_v5_1(z, k, gamma, false);`
 
Files
-----
`gsm_v5_1.m`      - Main Matlab function  
`gsm_v5_1_mex.c`  - Mex C code  
`README.md`       - This readme  
 
`gsm_v5_1_mex.mexw64`  - Mex binary for Windows  
`gsm_v5_1_mex.mexa64`  - Mex binary for Linux  
`makeit.m`             - Script to compile the Mex binary  
 
`./evaluation`         - Code for evaluating the numerical accuracy and running times.  
`./fast_math_binaries` - Mex binaries, compiled with the `fast-math` flag. Not used by default.
 
Evaluation of numerical accuracy and running times
--------------------------------------------------
CD to the `./evaluation` folder.  
  
To run the numerical accuracy test, assign some positive integer to the `task_id` variable,  
which determines the random seed, and call one of the test scripts. e.g.  
`>> task_id = 123;`  
`>> runTest2_double;`
 
The result is returned in a struct named `result`.  
Pre-run results over 200 instances are found in the corresponding `resultTest#_{double/single}.mat` files.  
 
To run the speed test, type  
`>> [times, dvals, kvals] = runSpeedTest.m;`
  
Pre-run results, averaged over 50 instances, are saved in `resultSpeedTest.mat`.
