# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 08:59:33 2022

@author: Sam Tuppen
"""

import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, freqz, firwin, lfilter
import segyio
import cvxpy as cp
import scipy.fftpack as spfft

#%% Parameters

dt = 0.004
flp = 30
fhi = 35
n = 7 #order Butterworth filter

tracenum = 99
end_time = 1024
start_time = 0


TracesPerGather = 635

filelocation = '/Users/samtu/Documents/'
filename = 'NoSI_short.segy'


#%% Import Data

with segyio.open(filelocation + filename, ignore_geometry=True) as f:
    #Load data
    k=0
    gather = f.trace.raw[:].T.reshape((len(f.samples),f.tracecount))[start_time:end_time,k+k*(TracesPerGather-1) : k+(k+1)*(TracesPerGather-1)]
    

#Add random noise
data = gather[0:end_time,tracenum] + 0.01*np.random.randn(*gather[0:end_time,tracenum].shape) # Select one trace and add noise
data = data - np.median(data)

#Defining some variables for FFT
Fs = 1/dt;
Nyq = Fs/2;
NFFT = int(2*len(data));  # Number of FFT
NFFT2 = int(NFFT/2)

left_l = np.ceil(20/Nyq*NFFT2)
left_r = np.ceil((flp-0.5)/Nyq*NFFT2)
right_l = np.ceil((fhi+0.5)/Nyq*NFFT2)
right_r = np.ceil(50/Nyq*NFFT2)
valid_index = np.hstack((np.arange(left_l,left_r,1), np.arange(right_l,right_r,1))).astype(int)

#%% Function
#Plot frequency and phase response
def mfreqz(b,a=1):
    w,h = freqz(b,a)
    h_dB = 20 * np.log10 (abs(h))
    plt.subplot(211)
    plt.plot(w/max(w),h_dB)
    plt.ylim(-150, 5)
    plt.ylabel('Magnitude (db)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Frequency response')
    plt.subplot(212)
    h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
    plt.plot(w/max(w),h_Phase)
    plt.ylabel('Phase (radians)')
    plt.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    plt.title(r'Phase response')
    plt.subplots_adjust(hspace=0.5)

#%% Filtering data

data_unfilt = data.copy()

Wn = [flp/Nyq, fhi/Nyq] #Butterworth non-dim freq
b, a = butter(n, Wn, 'bandstop', output='ba') #construct the filter
data = filtfilt(b, a, data, axis=0) #zero phase filter

#Generate a 'q', which is assumed to be the known wavelet
nq = 120
q = firwin(nq, cutoff = [2/250, 240/250], window = 'hamming', pass_zero = "bandpass")

# F, H_data = freqz(data, worN = NFFT2, fs = Fs)
# _, H_data_unfilt = freqz(data_unfilt, worN = NFFT2, fs = Fs)
# _, Hq = freqz(q, worN = NFFT2, fs = Fs) 

H_data = fft(data, n=NFFT, axis = 0)[0:NFFT2]
H_data_unfilt = fft(data_unfilt, n=NFFT, axis = 0)[0:NFFT2]
Hq = fft(q, n=NFFT, axis = 0)[0:NFFT2]


dftm = fft(np.eye(NFFT))
dftm = dftm[0:NFFT2, 0:NFFT]

div_factor = H_data[valid_index]/(Hq[valid_index]+0.00005)

#%% CVX
G0 = cp.Variable(NFFT)

objective = cp.Minimize(cp.norm(G0[0:NFFT2-1],1))
# HG = dftm@G0
constraints = [cp.norm((dftm@G0)[valid_index]-div_factor, 2) <= 0.00005]
prob = cp.Problem(objective, constraints)
prob.solve(verbose=True, solver = cp.ECOS, max_iters = 10)
# prob.solve(verbose=True, use_indirect=True, max_iters = 10000)


#%% Results

G = G0.value[0:NFFT2]

d_temp = lfilter(q, 1, G)
d = d_temp[0:NFFT2]

# _, HG = freqz(G, worN = NFFT2, fs = Fs)
# F, H_d = freqz(d, worN = NFFT2, fs = Fs)
HG = fft(G, n = NFFT2, axis = 0)
H_d = fft(d, n = NFFT2, axis = 0)
F = np.linspace(0,Nyq,NFFT2)

#%% Replacing the frequencies/phase to compute time
n=7
b, a = butter(n, Wn, 'bandpass', output='ba') #construct the filter
temp = filtfilt(b, a, d, axis=0) #zero phase filter
d_new = data + temp
H_new = fft(d_new, n=NFFT, axis=0)[0:NFFT2]
# H_new[int(left_r):int(right_l)] = H_d[int(left_r):int(right_l)]
# data_new = ifft(H_new)

#%% Plotting
plt.figure(figsize= (12,14))

plt.subplot(4,1,1)
plt.title('Time')
plt.plot(np.arange(0,dt*end_time-dt,dt),data_unfilt/np.std(data_unfilt),'m', label='Perfect')
plt.plot(np.arange(0,dt*end_time-dt,dt),d_new/np.std(d_new),'--k', label='Reconstruction')
plt.xlabel('Time')
plt.grid()
plt.xlim([0, dt*end_time-dt])
plt.legend()

plt.subplot(4,1,2)
plt.title('Amplitude')
plt.plot(F,20*np.log10(abs(H_data_unfilt)),'m', label='Perfect')
plt.plot(F,20*np.log10(abs(H_data)),'r', label='Data')
plt.plot(F,20*np.log10(abs(H_new)),'--b', label='Output')
# plt.plot(F,20*np.log10(abs(H_new)),'.k', label='Reconstruction')
plt.axvspan(left_l/NFFT2*Nyq, left_r/NFFT2*Nyq, color='green', alpha=0.1)
plt.axvspan(right_l/NFFT2*Nyq, right_r/NFFT2*Nyq, color='green', alpha=0.1)
plt.xlabel('Hz')
plt.xlim([20, 50])
plt.grid()
plt.legend()

plt.subplot(4,1,3)
plt.title('Phase')
plt.plot(F,np.angle(H_data_unfilt,deg=True),'m', label='Perfect')
plt.plot(F,np.angle(H_new,deg=True),'--b', label='Output')
plt.axvspan(left_l/NFFT2*Nyq, left_r/NFFT2*Nyq, color='green', alpha=0.1)
plt.axvspan(right_l/NFFT2*Nyq, right_r/NFFT2*Nyq, color='green', alpha=0.1)
plt.xlabel('Hz')
plt.ylabel('Degrees')
plt.xlim([20, 50])
plt.ylim([-180, 180])
plt.grid()
plt.legend()

plt.subplot(4,1,4)
plt.title('Phase error')
plt.stem(F,abs((np.angle(H_data_unfilt,deg=True))-(np.angle(H_new,deg=True))),'b', label='Error')
plt.axvspan(left_l/NFFT2*Nyq, left_r/NFFT2*Nyq, color='green', alpha=0.1)
plt.axvspan(right_l/NFFT2*Nyq, right_r/NFFT2*Nyq, color='green', alpha=0.1)
plt.axhline(y=10, color='r', linestyle='-')
plt.xlabel('Hz')
plt.ylabel('Degrees')
plt.xlim([20, 50])
plt.legend()
plt.ylim([0, 50])
plt.grid()
plt.show()



