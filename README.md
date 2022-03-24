# smbd_algorithms
Algorithms for Sparse Multichannel Blind Deconvolution (SMBD).

In this respository we share three different algorithms for the SMBD problem.

## c_pef.m

C-PEF algortihm

A cascade of a forward and a backward prediction error filter, optimized by means of a gradient descent algorithm to minimize the <img src="https://render.githubusercontent.com/render/math?math=\ell_1">-norm of the prediction error.

## am_smbd.m

AM-SMBD algorithm

An alternating minimization algorithm that alternates between the reflectivity series estimation (FISTA algorithm) and the wavelet estimation (Ridge regression estimator). Both algorithms are performed in the frequency domain.

## smbd_fista

FISTA-SMBD algorithm

The sparse multichannel blind deconvolution method (Euclid deconvolution + <img src="https://render.githubusercontent.com/render/math?math=\ell_1">-norm regularization)solved by the FISTA algorithm.

## main

Script that runs the algorithms for a synthetic dataset (synthetic2.mat).

## Authors

João Marcos Travassos Romano, PhD - UNICAMP 

Kenji Nose Filho, PhD - UFABC kenji.nose@ufabc.edu.br

Renan Del Buono Brotto, MSc - UNICAMP renanbrotto@gmail.com

Renato da Rocha Lopes, PhD - UNICAMP rlopes@unicamp.br

Thonia Cardoso Senna, MSc - UNICAMP 

## Documentation

See functions documentation and:

(To be published)
