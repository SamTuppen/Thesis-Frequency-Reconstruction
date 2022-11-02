%% Sparse solver
clear; %close all; clc;


%This script interpolates the missing energy in the time domain.

%% Data preparation

% Defining variables
timelimit = 1024;   % Remove last section of the trace (default should be about 8 for Pluto data)
dt = 0.004;         % Sampling time of the gather
dx = 12.5;          % Spatial sampling
c0 = 1400;          % must not be 0!
flp = 40;           % Lower bound of band-stop filter
fhi = 45;           % Upper bound of the band-stop filter

%Number of traces you with the interpolate for. As the script is slow for
%large amounts of data, it might be beneficial to compute only as many
%traces as are necessary.
tracelim = 400;     %The original gather has 636 traces.

%Selecting the dataset
file = "NoSI_200.segy";
filelocation = "C:\Users\samtu\OneDrive\Documents\University\Thesis\";

%Load data
[gather,STH,SH] = ReadSegy(filelocation + file);

%Add random noise
data = gather(1:timelimit,1:tracelim) + 0.01*randn(size(gather(1:timelimit,1:tracelim))); % Select one trace
data = data-median(data);

%Creating a new variable that will be bandstop filtered
datatemp = data;

%Define frequency param
Fs = 1/dt;              % Sampling frequency
Nyq = Fs/2;             % Nyquist frequency
NFFT = 2*length(datatemp);  % Number of FFT
NFFT2 = NFFT/2;

%% Filtering

datatemp = myFilter(datatemp,flp,fhi,dt);

rng(200) %Seed selected to add the same random noise. This corresponds to the noise added in the NN.
datatemp = datatemp + 0.1*randn(size(datatemp));


% Matching region for subject to condition in optimization, this is where
% the fit needs to match the known frequencies. So ideally, around the
% missing part. Hence, we have 2 sets below.
% They are shaped like so:
%
%      (known)    (unknown)    (known)
%   |------------|xxxxxxxxx|------------|
% left_l      left_r    right_l      right_r

left_l = ceil(20/Nyq*NFFT2); %frequency at which the solver needs to start matching the gather
left_r = ceil((flp-0.5)/Nyq*NFFT2); %Adding 0.5 Hz on each side to give the solver some 'breathing room'
right_l = ceil((fhi+0.5)/Nyq*NFFT2);
right_r = ceil(80/Nyq*NFFT2);   %frequency at which the solver needs to stop matching the gather

valid_index = [left_l:left_r, right_l:right_r];


%%
% Generate a 'q', which is assumed to be the known wavelet
q = fir1(120, [2 240]/250,'bandpass');              % Assumed time signature of the source wavelet (formula 1)
[Hq, ~] = freqz(q, 1, NFFT2, Fs);                   % Calculate source spectrum until Nyquist

% Create frequency transform matrix because CVX doesn't allow to do it later
dftm = dftmtx(NFFT);                % FFT matrix
dftm = dftm(1:NFFT2, 1:NFFT);       % FFT matrix until Nyquist (NFFT/2) for a vector of length NFFT

len_data=length(datatemp);          % Remember we selected NFFT2 as data length
d_recon = zeros(size(data));

for trace = 1:size(data,2)
    disp(trace)     %To track progress
    
    H_data = fft(datatemp(:,trace),NFFT);
    H_data = H_data(1:NFFT2);
    
    cvx_begin
        cvx_solver sedumi
        variables G0(NFFT,1);  % This is twice length of data
        minimize(norm(G0(1:NFFT2),1));   % minimize the L1 norm error (It should not be minimizing the L1 norm of the specter) (apply the norm on the first half of G0
        HG=dftm*G0; % Calculate FFT until Nyquist
        subject to
        norm(HG(valid_index)-H_data(valid_index)./(Hq(valid_index)+0.00005), 2) <= 0.00005;  % Between valid index frequencies
    cvx_end
    
    G=G0(1:len_data); %G is will be the reflectivity
    [HG, ~]=freqz(G,1,NFFT2,Fs);
    d=filter(q,1,G); %So then d is our approx. of the original signal
    d_cvx = d(1:len_data);
    
    %% Recon
    %Here we add the approximated energy to the frequency spectrum of the
    %filtered signal.

    d_temp = myFilter(d_cvx,flp,fhi,dt);

    d_filt = d_cvx - d_temp;
    
    d_recon(:,trace) = datatemp(:,trace) + d_filt;

end

%If you wish to save the data for later use and comparison then uncomment:
% save("gather200_Time_recon_40_50.mat",'d_recon','datatemp','data');




