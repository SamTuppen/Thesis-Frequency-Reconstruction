% Compressed Sensing
%
% Download and install CVX to run this example: http://cvxr.com/cvx/download/

% path(path,'C:\Users\kiranpour\Documents\SVD\Compressed_sensing_code\compressedsensing');
% path(path,'C:\Users\kiranpour\Documents\sam\Generate_gathers');
%
clear
close all;

% [Data,STH,SH]=ReadSegy();
% file = "SP_VQH_044_SI_bfSINAT_S1C12.segy";
file = "NoSI.segy";
filelocation = "C:\Users\samtu\Documents\";

%% Data preparation
trace = 120;
timelimit = 1024; %Remove last section of the trace (default should be about 8 for Pluto data)

%Load data
i=1;
[gather,STH,SH] = ReadSegy(filelocation + file,'traces',[i+(i-1)*635 : i+i*635]);


%Add random noise
data = gather(1:timelimit,trace) + 0.01*randn(size(gather(1:timelimit,trace))); % Select one trace
data = data-median(data);


%Define frequency param
dt = 0.004;
Fs=1/dt;
Nyq=Fs/2;
NFFT=2*length(data);  % Number of FFT
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
% left_l: left_r, right_l:right_r
left_l = ceil(20/Nyq*NFFT2);
left_r = ceil((flp-0.5)/Nyq*NFFT2);
right_l = ceil((fhi+0.5)/Nyq*NFFT2);
right_r = ceil(50/Nyq*NFFT2);
valid_index=[left_l:left_r, right_l:right_r];
% %Test
% test = round(1*length(valid_index));
% indices = randperm(length(valid_index));
% indices = indices(1:test);
% valid_index = valid_index(indices);


%% Filtering / Quadrature Mirror Filter (needs to be added)
n=7; %order Butterworth filter
Wn=[flp/Nyq fhi/Nyq]; %Butterworth non-dim freq
[b,a] = butter(n,Wn,'stop'); %construct the filter



datasaved = data;
data = filtfilt(b,a,data);   % Filter the trace, choose one wavelet
% data = StrongFilter(data,flp,fhi,dt);



%Generate a 'q', which is assumed to be the known wavelet
q=fir1(120, [2 240]/250,'bandpass');  % Assumed time signature of the source wavelet (formula 1)

[H_data, ~] = freqz(data, 1, NFFT2, Fs);  % Calculate data spectrum until Nyquist
[H_datasaved, ~] = freqz(datasaved, 1, NFFT2, Fs);  % Calculate data spectrum until Nyquist
[Hq, ~] = freqz(q, 1, NFFT2, Fs);  % Calculate source spectrum until Nyquist


%Create frequency transform matrix because CVX doesn't allow to do it later
dftm = dftmtx(NFFT);   % FFT matrix
dftm = dftm(1:NFFT2, 1:NFFT);  % FFT matrix until Nyquist (NFFT/2) for a vector of length NFFT

len_data=length(data);  % Remember we selected NFFT2 as data length


size(H_datasaved)
size(Hq)


%% Compressed Sensing
cvx_begin
%     cvx_solver ecos
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


%% Taking estimates and completing frequency spectrum

H_new = H_datasaved;
H_new(left_r:right_l) = H_d(left_r:right_l);

d_new = real(fft([H_new; flipud(H_new)]));
d_new = d_new(1:length(data));

%% Plotting
lw = 1;

%Figure showing what the trace looks like in the gather
%First plot shows where we are looking in the gather
figure();
subplot(4,2,1);
imagesc(gather)
colormap('gray')
caxis([-5 5])
hold on
xline(trace,'r','LineWidth',lw)
title('Trace being considered')

%Time Domain
subplot(4,2,2);
plot(datasaved/std(datasaved),'color','#77AC30','LineWidth',lw);  % Deconvolved using the above optimization
hold on;
% plot(d/std(d),'--b','LineWidth',lw);  % Deconvolved using the above optimization
plot(d_new/std(d_new),'r','LineWidth',lw)
grid on;
legend('Perfect','New','location','southeast');
xlabel('Samples');
xlim([250,350])
title('Time domain')

%Amplitude
subplot(4,2,[3,4]);
plot(F(1:end),20*log10(abs(H_d(1:end))),'--r','LineWidth',lw)
hold on
plot(F(1:end),20*log10(abs(H_datasaved(1:end))),'color','#77AC30','LineWidth',lw)
plot(F(1:end),20*log10(abs(H_data(1:end))),'k','LineWidth',lw);
grid on;
legend('Output','Data','Perfect','location','southeast');
xlabel('Hz');
ylim([-20 80])
%Adding highlights to plot
yl = ylim;
xl = xlim;
xBox = [left_l left_r left_r left_l;
        right_l right_r right_r right_l]'/NFFT2*Nyq;
yBox = [yl(2) yl(2) yl(1) yl(1);
        yl(2) yl(2) yl(1) yl(1)]';
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.03);
xlim([flp-10 fhi+10]);
set(gca,'children',flipud(get(gca,'children')),'YMinorTick','on')
title('Amplitude')

%Phase
subplot(4,2,[5,6]);
% plot(F(1:end),angle(H_d(1:end))*180/pi,'--b','LineWidth',lw)
plot(F(1:end),angle(H_new(1:end))*180/pi,'--r','LineWidth',lw)
hold on
% plot(F(1:end),angle(H_data(1:end))*180/pi,'r','LineWidth',lw);
plot(F(1:end),angle(H_datasaved(1:end))*180/pi,'color','#77AC30','LineWidth',lw)
grid on;
% legend('Output','Data','Perfect','location','southeast');
xlabel('Hz');
ylim([-180 180])
%Adding highlights to plot
yl = ylim;
xl = xlim;
xBox = [left_l left_r left_r left_l;
        right_l right_r right_r right_l]'/NFFT2*Nyq;
yBox = [yl(2) yl(2) yl(1) yl(1);
        yl(2) yl(2) yl(1) yl(1)]';
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.03);
xlim([flp-10 fhi+10]);
set(gca,'children',flipud(get(gca,'children')),'YMinorTick','on')
legend('Provided data','Perfect','New')
title('Phase')

%Phase error
subplot(4,2,[7,8]);
stem(F(1:end),abs(angle(H_datasaved(1:end))-angle(H_d(1:end)))*180/pi,'Color',[0 0.498039215803146 0],'LineWidth',lw);
hold on;
stem(F(1:end),abs(angle(H_datasaved(1:end))-angle(H_new(1:end)))*180/pi,'Color','red','LineWidth',lw);
grid on;
% legend('Output','Data','Perfect','location','southeast');
xlabel('Hz');
%Adding highlights to plot
yl = ylim;
xl = xlim;
xBox = [left_l left_r left_r left_l;
        right_l right_r right_r right_l]'/NFFT2*Nyq;
yBox = [yl(2) yl(2) yl(1) yl(1);
        yl(2) yl(2) yl(1) yl(1)]';
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.03);
yline(10,'LineWidth',lw)
xlim([flp-10 fhi+10]);
set(gca,'children',flipud(get(gca,'children')),'YMinorTick','on')
title('Phase error')

