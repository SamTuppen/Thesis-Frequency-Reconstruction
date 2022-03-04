function [gather, gatherfilt] = create_linear_synth(Number_of_traces, Number_of_samples, dt, dx, v, trace_offset, sample_offset)
% gather=create_linear_synth(Number_of_traces, Number_of_samples, v, trace_offset, sample_offset)
%  
%  dt = Sampling period (2ms)
%  dx = Distance between sensors in meters also called spatial sampling (m)
%  v is velocity
%  trace_offset is some offset in number of traces where signal starts 
%  sample_offset is some offset to samples where signal starts

%% Wavelet parameters
Fs=1/dt;    % Sampling frequency (500 Hz)
Nyq=Fs/2;   % Nyquist frequency
NFFT=8192;  % Number of fft

[b,a]=fir1(15,0.8);  % Create an FIR filter (a returns always 1 for FIR filters)
b=b(5:end)-0.99*b(1:end-4); % Create ghost

%% Plot wavelet
% figure();clf;
% plot(b);grid on;
% ylabel('Amplitude');
% xlabel('Coefficient number')
% title('My wavelet');
% 
% figure();clf;
% freqz(b,a,NFFT,Fs);  % Frequency response of your filter
% title('Ghost''s frequency response');

%% Generate line

% Use b as a wavelet for your synthetic data
length_b=length(b);

gather=zeros(Number_of_samples,Number_of_traces); % Create gather


for ii=trace_offset+1:Number_of_traces  % Create a gather with two moveouts
    from_offset=ii-trace_offset;
    entire_sample=round(from_offset*dx/v*Fs);
    part = from_offset * dx / v * Fs - entire_sample;
    entire_sample = entire_sample+sample_offset;
    gather(entire_sample : entire_sample - 1 + length_b,ii) = b;  % Move entire samples
    gather(:,ii) = partial_sample_shift(gather(:,ii),part);  % Move partial
end

% figure(3);
% imagesc([0 Number_of_traces-1]*dx,[0 Number_of_samples-1]*dt,gather);
% xlabel('Meter'); 
% ylabel('Seconds');
% title('My simple synthetic');
% %caxis([-0.2 0.2]);
% colorbar;
% colormap('gray'); % Or some other colormap

gather = gather(1 : Number_of_samples , 1 : Number_of_traces);


