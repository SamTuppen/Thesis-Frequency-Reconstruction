# README - Frequency Band Interpolation

Two different methods to interpolate missing frequency bands are presented in this folder. The first is an approach using a neural network (NN) and the second is a sparse solver (SS) approach. The NN has been created in Python using the TensorFlow platform and the SS has been created in MATLAB.  

---

In their current state, the following methods are set to interpolate the frequency data of a dataset called `NoSI_short.segy`, which is a cropped version of the original `SP_VQH_044_S1C12.segy` dataset.
Furthermore,`CroppingInputData.m` has been included to generate the training and test samples for the neural network. It will be necessary to run this script prior to running the network.

# NEURAL NETWORK SECTION:

To run the NN you will need (Python):
- TensorFlow
- NumPy
- Matplotlib
- SciPy
- tqdm (for plotting the loading bars)
- segyio (for loading the seismic data)

_(To run the neural network on your graphics card you will need to install additional software, for example CUDA, that is available to find in the TensorFlow documentation: https://www.tensorflow.org/install/pip)_

Several minor variations of the script exist, but the main code has been included: `Neural Network Interpolation.ipynb`.

Structure of the code:
- Import packages
- Define necessary parameters
- Define functions
- Load data and generate training/test samples
- Build Neural Network
- Fit the network to training data
- Make predictions
- QC

Within the above script, a CA-Unet section can be found, which is a Coordinate Attentive block that replaces the encoding blocks in the standard U-Net. These blocks increase the training time, but help the network become more spatially aware of surrounding traces. In general I did not observe a great improvement, but I have left them in the code (commented out) in case anyone is interested in trying them. Original code in pyTorch: https://github.com/Andrew-Qibin/CoordAttention.  
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

To run the MATLAB code you will need to install CVX, which can be downloaded using the following link: 
http://cvxr.com/cvx/download/  
You will also need [SegyMAT](http://segymat.sourceforge.net/), which is used to read the seismic data in MATLAB.

There are several scripts:
1. `Sparse_Solver_Time_Domain.m`: Reconstruct purely in the time domain
2. `Sparse_Solver_Taup_Domain.m`: Reconstructs in the tau-p domain
3. `gather200_timeRecon.mat`: This is a data file of a time domain reconstruction with 10% noise
4. `gather200_timeRecon_no_noise.mat`: This is a data file of a time domain reconstruction
5. `gather200_TaupRecon.mat`: This is a data file of a taup domain reconstruction with 10% noise
6. `taup_example.m`: Example file of taup transform
7. `NoSI_200.segy`: 200th shot gather of the above metioned dataset.

The MATLAB interpolation scripts are structured as follows:
1. Define survey parameters
2. Specify file location
3. Load the data
4. Specify which frequencies are known and which are unknown
5. Filter the data to create missing frequency band that must be interpolated
6. Compute sparse 'reflectivity' of the gather
7. The interpolated band is added to the sparse spectrum to complete the spectrum


