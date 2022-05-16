# README - Frequency Interpolation (Sam Tuppen)

There are two different sections for interpolating missing frequency bands. The first is an approach using a neural network (NN) and the second is the sparse solver (SS) approach. The NN has been created in Python using the TensorFlow platform and the SS has been created in both Python and MATLAB.

# NEURAL NETWORK SECTION:

To run the NN you will need (Python):
- TensorFlow
- NumPy
- Matplotlib
- SciPy
- tqdm (for plotting the loading bars)

_(To run the network on your graphics card you will need to install additional software, for example CUDA, that is available to find in the TensorFlow documentation.)_


# SPARSE SOLVER SECTION:

## MATLAB:
To run the MATLAB code you will need to install CVX, which can be downloaded using the following link: 
http://cvxr.com/cvx/download/

There are several scripts:
1. `TVNorm_SparseReg.m`: Reconstruct one missing frequency band for one trace
2. `TVNorm_SparseReg_multiple_gaps.m`: Reconstructs multiple missing frequency bands for one trace


## Python:
To Run the Python code you will need the following packages:
- cvxpy
- NumPy
- Matplotlib
- SciPy
