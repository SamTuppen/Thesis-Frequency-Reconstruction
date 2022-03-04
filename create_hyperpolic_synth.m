function [gather, gatherfilt] = create_hyperpolic_synth(Number_of_traces, Number_of_samples, dt, dx, v, h, xshift)
% gather=create_linear_synth(Number_of_traces, Number_of_samples, v, trace_offset, sample_offset)
%  
%  v is velocity
%  h is depth of reflector in metres
%  trace_offset is some offset in number of traces where signal starts 
%  sample_offset is some offset to samples where signal starts
%  dt = Sampling period (2ms)
%  dx = Distance between sensors in meters also called spatial sampling (m)



% Number_of_traces = 128;
% Number_of_samples  = 256;
% dt = 0.002;
% dx = 6.25;
% v = 1.303666137722284e+03;
% h = 10;




%% Parameters
Fs=1/dt;    % Sampling frequency (500 Hz)
Nyq=Fs/2;   % Nyquist frequency
NFFT=8192;  % Number of fft


%% Create wavelet

[b,a]=fir1(15,0.8);         % Create an FIR filter (a returns always 1 for FIR filters)
b=b(5:end)-0.99*b(1:end-4); % Create ghost

% figure();clf;
% plot(b);grid on;
% ylabel('Amplitude');
% xlabel('Coefficient number')
% title('My wavelet');
% 
% figure();clf;
% freqz(b,a,NFFT,Fs);  % Frequency response of your filter
% title('Ghost''s frequency response');


%% Generate gather

% Use b as a wavelet for your synthetic data
length_b=length(b);

%Add buffers to make smooth edges
Number_of_traces = Number_of_traces + length_b;
Number_of_samples = Number_of_samples + length_b;

%Create distance and time axis
x = [0:dx:Number_of_traces * dx - dx] - xshift;
t = [0:dt:Number_of_samples * dt - dt];

%create matrix that will be filled in
gather=zeros(Number_of_samples,Number_of_traces); % Create gather

%% Initial attempt

% %Fill in matrix based on physics
% for ii = 1:Number_of_traces
%     
%     tIndex = sqrt(x(ii).^2 + 4*h^2)/v;
%     [~,it] = min(abs(t - tIndex));
% 
%     if it+ceil(length_b/2) > size(gather,1)
%         break
%     end
%     gather(it-floor(length_b/2) + 1 : it+ceil(length_b/2) , ii) = b;
% end

%% Second attempt (Kambiz Linear event code)

% for ii = xshift+1 : Number_of_traces  % Create a gather with two moveouts
%     from_offset = ii - xshift;
%     entire_sample = round(from_offset*dx/v*Fs);
%     part = from_offset * dx / v * Fs - entire_sample;
%     entire_sample = entire_sample + sample_offset;
% 
%     gather(entire_sample : entire_sample - 1 + length_b,ii) = b;  % Move entire samples
%     gather(:,ii) = partial_sample_shift(gather(:,ii),part);  % Move partial
% end


%% Third attempt

%Fill in matrix based on physics
for ii = xshift+1 : Number_of_traces


    from_offset = ii - xshift;
    entire_sample = round(from_offset*dx/v*Fs);
    part = from_offset * dx / v * Fs - entire_sample;

    tIndex = sqrt(x(ii).^2 + 4*h^2)/v;
    [~,it] = min(abs(t - tIndex));

    if it+ceil(length_b/2) > size(gather,1)
        break
    end
    
    gather(it-floor(length_b/2) + 1 : it+ceil(length_b/2) , ii) = b;
    gather(:,ii) = partial_sample_shift(gather(:,ii),part);  % Move partial
end


%Remove buffer
gather = gather(1 : end-length_b , 1 : end-length_b);
