%%Main File
clear; close all; clc;

%% Defining parameters

Ntraces = 128;              % Number of traces [-]
Nsamples = 256;             % Number of time measurements [-]
traceOffset = 6.25;         % Distance between traces [m]
Ntrainingsamples = 100;       % Number of training examples to generate [-]
dt = 0.002;                 % Sampling rate [s]
dx = 6.25;                  % Spatial Sampling distance [m]
NFFT = 8192;                % Number of fft
Fs = 1/dt;                  % Sampling frequency (500 Hz)
df = 0.1;                   % Frequency resolution in Hz
xshift = 0;                 %Amount to shift hyperbolas in x-direction [m]
noiselevel = 0.1;           %Percent noise
freqRemove = [];            %Frequency range to be removed [Hz]; Parameter ignored in case of random ranges
RandomFilts = 1;            %Turn on for using 2 random 10 Hz frequency ranges

SaveResults = 1;           %Save plots? 1: Yes,
%                                        0: No.
SaveNormFac = 1;            %Save Normalization factors? Needed for fitting data later

% FileLocationX = 'C:\Users\samtu\Documents\XDataTrain';         %Where to save Training data X
% FileLocationY = 'C:\Users\samtu\Documents\YDataTrain';         %Where to save Training data Y
FileLocationX = 'C:\Users\samtu\Documents\XDataTest';         %Where to save Training data X
FileLocationY = 'C:\Users\samtu\Documents\YDataTest';         %Where to save Training data Y

%% Generate gather with hyperbolic and linear events
fprintf('Creating synthetic gathers...\n')
for i = 1:Ntrainingsamples


    %Add linear event
    v = 1300 + (4000 - 1300) .* rand(1,1);
    trace_offset = 0;                                    %Trace offset to start drawing linear event
    sample_offset = ceil(0 + (Nsamples-10 - 0) .* rand(1,1));          %Time offset to start drawing linear event
    gatherTemp1 = create_linear_synth(Ntraces, Nsamples, dt, dx, v, trace_offset, sample_offset);

    %Add hyperbolic event 1
    v = 1300 + (4000 - 1300) .* rand(1,1);
    h = 20 + (450 - 20) .* rand(1,1);
    gatherTemp2 = create_hyperpolic_synth(Ntraces, Nsamples, dt, dx, v, h, xshift);

    %Add hyperbolic event 2
    v = 1300 + (4000 - 1300) .* rand(1,1);
    h = 20 + (450 - 20) .* rand(1,1);
    gatherTemp3 = create_hyperpolic_synth(Ntraces, Nsamples, dt, dx, v, h, xshift);

    %Put events together into one matrix
    gather = gatherTemp1 + gatherTemp2 + gatherTemp3;

    %Progress counterimage
    if mod(i,100) == 0
        fprintf('Progress: %d %% \n', round(i/Ntrainingsamples*100))
    end

    %% Generate Filtered Gather

    %Filter
    if RandomFilts == 1
    freqRemove(i,1) = round(1 + (Fs/2-20 - 1) .* rand(1,1)); %Filtering out the a random frequency range with a bandwidth of 10 Hz. I left the last 10 Hz untouched by the filters as this range simply includes edge effects that are not of interest. Interpolation will become too difficult and it it not a realistic scenario
    freqRemove(i,2) = freqRemove(i,1) + 10;

    gatherfilt = myFilter(gather,freqRemove(i,1),freqRemove(i,2),dt);

    freqRemove(i,3) = round(1 + (Fs/2-20 - 1) .* rand(1,1)); 
    freqRemove(i,4) = freqRemove(i,3) + 10;

    gatherfilt = StrongFilter(gatherfilt,freqRemove(i,3),freqRemove(i,4),dt); %Removing 2 random bands

    else
    gatherfilt = StrongFilter(gather,freqRemove(1),freqRemove(2),dt); %Can only remove one frequency band when picking a range

    end
%     PlotSpectra(gather, gatherfilt, Ntraces, Nsamples, dx, dt)

    %Add Noise
    gatherfilt = gatherfilt + noiselevel*randn(size(gatherfilt));
   

    %% Normalize to [0 1]
    gatherfactor(i) = max(max(abs(gather)));
    gatherfiltfactor(i) =  max(max(abs(gatherfilt)));

    gather = ((gather ./ max(max(abs(gather)))) + 1) ./ 2;

    gatherfilt = ((gatherfilt ./ max(max(abs(gatherfilt)))) + 1) ./ 2;


    %Check to see if lowest value is not less than 0.
    if min(min(gather)) < 0 || min(min(gatherfilt)) < 0
        error('LESS THAN ZERO! Current iteration: %d',i)
    end

    %% Plot Gathers if necessary

%     PlotGathers(gather, gatherfilt, Ntraces, Nsamples, dx, dt);

    %% Save results

    if SaveResults == 1

        %Save unfiltered (X Data)
        FileNameX = FileLocationX + "\" + "X_" + string(i);
        writematrix(gatherfilt,FileNameX)
        %         imwrite(gatherfilt, FileNameX);
        %         save(FileNameX,'gather')

        %Save filtered (Y Data)
        FileNameY = FileLocationY + "\" + "Y_" + string(i);
        writematrix(gather,FileNameY)
        %         imwrite(gather, FileNameY);
        %         save(FileNameY,'gatherfilt')

    end
   
end

if SaveNormFac == 1

    writematrix(gatherfactor,'gatherfactor');
    writematrix(gatherfiltfactor,'gatherfiltfactor');
    
    writematrix(freqRemove,'freqRemove');
end

fprintf('Done!\n')
