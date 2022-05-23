# README - Frequency Interpolation (Sam Tuppen)

There are two different sections for interpolating missing frequency bands. The first is an approach using a neural network (NN) and the second is the sparse solver (SS) approach. The NN has been created in Python using the TensorFlow platform and the SS has been created in both Python and MATLAB.

# NEURAL NETWORK SECTION:

To run the NN you will need (Python):
- TensorFlow
- NumPy
- Matplotlib
- SciPy
- tqdm (for plotting the loading bars)

_(To run the neural network on your graphics card you will need to install additional software, for example CUDA, that is available to find in the TensorFlow documentation: https://www.tensorflow.org/install/pip)_

There are several scripts:
1. `Full U-Net FBW.ipynb`: NN code with Filtering Before Windowing (FBW)
2. `Full U-Net FAW.ipynb`: NN code with Filtering After Windowing (FAW), this code is mainly used when we wish to completely zero the missing frequency bands.
3. 

Note: I intend to make these into `.py` files once everything runs correctly.  
Several of the above scripts have a CA-Unet sections, with is a Coordinate Attentive block that replaces the encoding blocks in the standard U-Net. These blocks increase the training time, but help the network become more spatially aware of surrounding traces. In general I did not observe a great improvement, but I have left them in the code (commented out) in case anyone is interested in trying them.

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


