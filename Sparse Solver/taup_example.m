%Testing the tau-p algorithm (in 2d) from Robertsson (Geophysics 2019).
%Finding: On a very fine grid - it is almost perfect
%On normal data - is is crappy --> need to interpolate a lot....

clear; close all;clc;


file = "NoSI.segy";
filelocation = "C:\Users\samtu\Documents\";
timelimit = 1500; %Remove last section of the trace (default should be about 8 for Pluto data)
i=1; %gather number
[gather,~,~] = ReadSegy(filelocation + file,'traces',[i+(i-1)*635 : i+i*635]);
data = gather(1:timelimit,:); % Select one trace

% [samps_ts, no_tr] = size(data);

a = 100;
f = data;

%n=100
%f=f(1+n:550+n,1:550);
% f=f(1:1500,:);
figure
subplot(1,4,1)
imagesc(f, [-a, a]); colormap gray
title(sprintf('Input: %g',sum(sum(abs(f)))))

%forward transform
dt=0.004;
dx=12.5;
c0=1000; %must not be 0!
[g,tau,p]=taup2(f,dt,dx,c0);

subplot(1,4,4)
imagesc(g, [-1000 1000]); colormap gray
title('taup shot')

%backward transform
[data_out,tau,p]=taup2inv(g,dt,dx,c0);
data_out=fliplr(data_out);

subplot(1,4,2)
imagesc(data_out, [-a a]); colormap gray
title(sprintf('Out: %g',sum(sum(abs(data_out)))))

subplot(1,4,3)
imagesc(f-data_out, [-a a]); colormap gray
title(sprintf('diff: %g',sum(sum(abs(f-data_out)))))

% %%Now a solution would be to interpolate to a much finer grid - in the
% %%x-dir 
% F=interp2(f,3,'cubic');
% %F=F(1:5990,1:2250);
% F=F(1:11990,1:4500);
% size(f)
% size(F)
% figure; 
% subplot(1,4,1)
% imagesc(F, [-a a])
% title(sprintf('Input: %g',sum(sum(abs(F)))))
% 
% colormap gray
% 
% %Forward transform
% dx=12.5/8; dt=0.004/8;
% [G,tau,p]=taup2(F,dt,dx,c0);
% 
% %backward transform
% [data_out,tau,p]=taup2inv(G,dt,dx,c0);
% data_out=fliplr(data_out);
% 
% subplot(1,4,4)
% imagesc(g, [-1000 1000]); colormap gray
% title('taup shot')
% 
% subplot(1,4,2)
% imagesc(data_out, [-a a]); colormap gray
% title(sprintf('Out: %g',sum(sum(abs(data_out)))))
% 
% subplot(1,4,3)
% imagesc(F-data_out, [-a a]); colormap gray
% title(sprintf('diff: %g',sum(sum(abs(F-data_out)))))







