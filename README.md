# README - Frequency Band Interpolation

Two different methods to interpolate missing frequency bands are presented in this folder. The first is an approach using a neural network (NN) and the second is a sparse solver (SS) approach. The NN has been created in Python using the TensorFlow platform and the SS has been created in both Python and MATLAB.  

Both methods have been tested for the 
---

In their current state, the following methods are set to interpolate the frequency data of a dataset called `NoSI_short.segy`, which is a cropped version of the original `SP_VQH_044_S1C12.segy` dataset.

# NEURAL NETWORK SECTION:

To run the NN you will need (Python):
- TensorFlow
- NumPy
- Matplotlib
- SciPy
- tqdm (for plotting the loading bars)
- segyio (for loading the seismic data)

_(To run the neural network on your graphics card you will need to install additional software, for example CUDA, that is available to find in the TensorFlow documentation: https://www.tensorflow.org/install/pip)_

There are several scripts:
1. `Full U-Net FBW.ipynb`: NN code with Filtering Before Windowing (FBW)
2. `Full U-Net FAW.ipynb`: NN code with Filtering After Windowing (FAW), this code is mainly used when we wish to completely zero the missing frequency bands.
3. 

Note: I intend to make these into `.py` files once everything runs correctly.  
Several of the above scripts have a CA-Unet sections, with is a Coordinate Attentive block that replaces the encoding blocks in the standard U-Net. These blocks increase the training time, but help the network become more spatially aware of surrounding traces. In general I did not observe a great improvement, but I have left them in the code (commented out) in case anyone is interested in trying them. Original code in pyTorch: https://github.com/Andrew-Qibin/CoordAttention.  
Furthermore, the paper in which I found them being used on seismic data is: [Li, X., Wu, B., Zhu, X., & Yang, H. (2021). Consecutively Missing Seismic Data Interpolation based on Coordinate Attention Unet. _IEEE Geoscience and Remote Sensing Letters_, 19, 1-5](https://ieeexplore.ieee.org/document/9615194).

The structure of the U-Net is as follows:  
<p align="center">
  <img 
    width="430"
    height="350"
    src="https://user-images.githubusercontent.com/93287046/169788922-bc895f53-c690-4836-9a4d-a4ad10a5b287.png"
  >
</p>


# SPARSE SOLVER SECTION:
The sparse solver approach is an implementation of the scheme presented by [Wang, R. & Herrmann, F.](https://doi.org/10.1190/segam2016-13879674.1).  
## MATLAB:
To run the MATLAB code you will need to install CVX, which can be downloaded using the following link: 
http://cvxr.com/cvx/download/  
You will also need [SegyMAT](http://segymat.sourceforge.net/), which is used to read the seismic data in MATLAB.

There are several scripts:
1. `TVNorm_SparseReg.m`: Reconstruct one missing frequency band for one trace
2. `TVNorm_SparseReg_multiple_gaps.m`: Reconstructs multiple missing frequency bands for one trace

The MATLAB scripts are structured as follows:
1. Define survey parameters
2. Specify file location
3. Load the data
4. Specify which frequencies are known and which are unknown
5. Filter the data to create missing frequency band that must be interpolated
6. Compute sparse 'reflectivity' of the gather
7. The interpolated band is added to the sparse spectrum to complete the spectrum


## Python:
To Run the Python code you will need the following packages:
- cvxpy
- NumPy
- Matplotlib
- SciPy


