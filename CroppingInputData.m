% clear;

%This is a simple script to crop the input data so that it is easier to
%load into the neural network

file = 'SP_VQH_044_S1C12.segy';

%% Load Data
%You can either load the full file, or only certain traces.

% [Data,STH,SH]=ReadSegy(file);


% In the following manner you can select different traces/gather:
% For example, to select the first 200 gathers, where each gather contrains 636
% traces use the following command:

[Data,STH,SH]=ReadSegy(file,'traces',[1:200*636]); 


WriteSegyStructure('NoSI_short.segy',SH,STH,Data); %Write data to a new fill
clear 
[Data,STH,SH]=ReadSegy('NoSI_short.segy');%Open the new file to test

%% Plotting
figure
imagesc(Data);
colormap('gray');
caxis([-1 1]);