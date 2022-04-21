% Compressed Sensing
%
% Download and install CVX to run this example: http://cvxr.com/cvx/download/

% path(path,'C:\Users\kiranpour\Documents\SVD\Compressed_sensing_code\compressedsensing');
% path(path,'C:\Users\kiranpour\Documents\sam\Generate_gathers');
%
% clear all
% close all;

%% Data preparation
GenerateTrainingFigures;

%Add random noise
data=gather(:,100)+0.01*randn(size(gather(:,100))); % Select one trace
data=data-median(data);

%Define frequency param
Fs=500;
Nyq=Fs/2;
NFFT=512;  % Number of FFT
NFFT2=NFFT/2;
NFFT4=NFFT/4;
NFFT6=NFFT/6;
NFFT8=NFFT/8;

%Missing frequencies
flp = 30;
fhi = 35;


% Matching region for subject to condition in optimization, this is where
% the fit needs to match the known frequencies. So ideally, around the
% missing part. Hence, we have 2 sets below.
%left_l: left_r, right_l:right_r
left_l = 1;
left_r = ceil(flp/Nyq*NFFT2);
right_l = ceil(fhi/Nyq*NFFT2);
right_r = ceil(NFFT4);
valid_index=[left_l: left_r, right_l:right_r];



%% Filtering / Quadrature Mirror Filter (needs to be added)
n=7; %order Butterworth filter
Wn=[flp/Nyq fhi/Nyq]; %Butterworth non-dim freq
[b,a] = butter(n,Wn,'stop'); %construct the filter

datasaved = data;
data = filtfilt(b,a,data(1:NFFT2));   % Filter the trace, choose one wavelet


%Generate a 'q', which is assumed to be the known wavelet
q=fir1(70, [2 240]/250,'bandpass');  % Assumed time signature of the source wavelet (formula 1)


[H_data, ~] = freqz(data, 1, NFFT2, Fs);  % Calculate data spectrum until Nyquist
[H_datasaved, ~] = freqz(datasaved, 1, NFFT2, Fs);  % Calculate data spectrum until Nyquist
[Hq, ~] = freqz(q, 1, NFFT2, Fs);  % Calculate source spectrum until Nyquist


%Create frequency transform matrix because CVX doesn't allow to do it later
dftm = dftmtx(NFFT);   % FFT matrix
dftm = dftm(1:NFFT2, 1:NFFT);  % FFT matrix until Nyquist (NFFT/2) for a vector of length NFFT

len_data=length(data);  % Remember we selected NFFT2 as data length


%% Compressed Sensing
cvx_begin
    variables G0(NFFT,1);  % This is twice length of data
    minimize(norm(G0(1:NFFT2),1));   % minimize the L1 norm error (It should not be minimizing the L1 norm of the specter) (apply the norm on the first half of G0
    HG=dftm*G0; % Calculate FFT until Nyquist
    subject to
    norm(HG(valid_index)-H_data(valid_index)./(Hq(valid_index)+0.00005), 2) <= 0.00005;  % Between valid index frequencies
cvx_end


G=G0(1:len_data);
[HG, ~]=freqz(G,1,NFFT/2,Fs);
d=filter(q,1,G);
d=d(1:len_data);

[H_d, F] = freqz(d,1,NFFT/2,Fs);


%% Plotting
lw = 1;

figure(1);
title('Time domain')
plot(datasaved/std(datasaved),'m','LineWidth',lw);  % Deconvolved using the above optimization
hold on;
plot(data/std(data),'r','LineWidth',lw);  % Original data
plot(d/std(d),'--b','LineWidth',lw);  % Deconvolved using the above optimization
grid on;
legend('Perfect','Data','Output','location','southeast');
xlabel('Samples');


figure(2);
title('Amplitude')
plot(F(1:end),20*log10(abs(H_d(1:end))),'--b','LineWidth',lw)
hold on
plot(F(1:end),20*log10(abs(H_data(1:end))),'r','LineWidth',lw);
plot(F(1:end),20*log10(abs(H_datasaved(1:end))),'m','LineWidth',lw)
grid on;
legend('Output','Data','Perfect','location','southeast');
xlabel('Hz');
ylim([-180 180])
%Adding highlights to plot
yl = ylim;
xl = xlim;
xBox = [flp fhi fhi flp];
yBox = [yl(2) yl(2) yl(1) yl(1);];
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.05);
xlim([flp-20 fhi+20]);
set(gca,'children',flipud(get(gca,'children')),'YMinorTick','on')

figure(3);
title('Phase')
plot(F(1:end),angle(H_d(1:end))*180/pi,'--b','LineWidth',lw)
hold on
plot(F(1:end),angle(H_data(1:end))*180/pi,'r','LineWidth',lw);
plot(F(1:end),angle(H_datasaved(1:end))*180/pi,'m','LineWidth',lw)
% stem(F(1:end),abs(angle(H_datasaved(1:end))-angle(H_d(1:end)))*180/pi,'LineWidth',lw);
grid on;
legend('Output','Data','Perfect','location','southeast');
xlabel('Hz');
ylim([-180 180])
%Adding highlights to plot
yl = ylim;
xl = xlim;
xBox = [flp fhi fhi flp];
yBox = [yl(2) yl(2) yl(1) yl(1);];
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.05);
xlim([flp-20 fhi+20]);
set(gca,'children',flipud(get(gca,'children')),'YMinorTick','on')
