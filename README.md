Code to calculate the Generalized Soft-Min function
===================================================

Version 5.1, June 2020

Based on the paper: 
Tal Amir, Ronen Basri, Boaz Nadler (2020) - The Trimmed Lasso: Sparse Recovery Guarantees and Practical Optimization by the Generalized Soft-Min Penalty

Requirements
------------
**Matlab**  
This program supports Matlab 2018b and onward, but may work with earlier versions.  
The .mex binaries are compiled for Windows and Linux. 
On Mac, the code can be compiled on the target computer, or run without the mex code; see below.
  
Usage
-----
Basic usage:
`>> [mu, theta] = gsm_v5_1(z, k, gamma);`

To avoid calling the Mex binary and run purely on Matlab:
`>> [mu, theta] = gsm_v5_1(z, k, gamma, false);`

Files
-----
`sparse_approx_gsm_v1_21.m`    - Main Matlab function  
`sparse_approx_gsm_v1_21.txt`  - Main documentation  
`README.md`                    - This readme  

`runExample*.m`             - Script files with simple usage examples  
`runCompareTrimmedLasso.m`  - A comparison between GSM and the DC-Programming and ADMM methods described in [2].
                          
`./gsm`  - Required in the matlab path.  
`./comparison`  - Required only for comparing GSM with other methods. 

Evaluation of numerical accuracy and running times
--------------------------------------------------
The `./evaluation` folder contains code that runs the numerical accuracy evaluation
and measures the running times.

To run the numerical accuracy test, assign some positive integer to the `task_id` variable,
which determines the random seed, and call one of the test scripts. e.g.
`>> task_id = 123;`
`>> runTest2_double;`

The result is returned in a struct named `result`.
Pre-run results over 200 instances are found in the corresponding `resultTest#_{double/single}.mat` files.

To run the speed test, type
`>> [times, dvals, kvals] = runSpeedTest.m;`

Pre-run results, averaged over 50 instances, are saved in `resultSpeedTest.mat`.
