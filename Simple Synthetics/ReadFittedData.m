clear; close all; clc;
%%

PlotFreq = 0;       %Plot the frequency and phase response
PlotTime = 0;       %Plot the data in the time-domain
CalcSens = 1;       %Calculate the 'Sensitivity'. (RMSE)


for index = 1:100

    Ntraces = 128;              % Number of traces [-]
    Nsamples = 256;             % Number of time measurements [-]
    dt = 0.002;                 % Sampling rate [s]
    dx = 6.25;                  % Spatial Sampling distance [m]
    Fs = 1/dt;
    NFFT = 8192;


    FileLocationX = "Neural Network\XDataTest\";
    FileLocationY = "Neural Network\YDataTest\";
    FileLocationZ = "Neural Network\ZDataTest\";

    gatherfactor = readmatrix("gatherfactor30" + ".txt");
    gatherfiltfactor = readmatrix("gatherfiltfactor30" + ".txt");

    DataX = readmatrix(FileLocationX + "X_" + string(index) + ".txt");
    DataY = readmatrix(FileLocationY + "Y_" + string(index) + ".txt");
    DataZ = readmatrix(FileLocationZ + "Z_" + string(index) + ".txt");
    
    fprintf('Interation Complete: %.f \n', index)

    %% Frequency Plotting
        input=DataY(:,1);
        filtered = DataX(:,1);
        out=DataZ(:,1);

        f1=abs(fft(input));
        f2=abs(fft(out));
        f3=abs(fft(filtered));
        x=linspace(0,500,length(out));
        
        p1 = angle(fft(input));
        p2 = angle(fft(out));

     if PlotFreq == 1

        %Load Filter boundaries
        freqRemove = readmatrix("freqRemove.txt");
        set(0,'defaulttextinterpreter','latex')

        figure();
        subplot(2,1,1)
        hold on

        plot(x,20*log10(f1),'r')
        plot(x,20*log10(f2),'b')
        plot(x,20*log10(f3),'--','color','#00841a')
        axis([0,250,-30,10])
        title('\textbf{Frequency response}');
        xlabel('\textbf{Frequency [Hz]}')
        ylabel('\textbf{Magnitude [dB]}')
        set(gca,'fontname','times')
        set(gca,'Fontsize',12);

        %Adding highlights to plot
        yl = ylim;
        xl = xlim;

        xBox = [freqRemove(index,1) freqRemove(index,3);
                freqRemove(index,1) freqRemove(index,3);
                freqRemove(index,2) freqRemove(index,4);
                freqRemove(index,2) freqRemove(index,4)];

        yBox = [yl(2) yl(2);
                yl(1) yl(1);
                yl(1) yl(1);
                yl(2) yl(2)];
        patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.1);

        grid on
        legend('Ground Truth','Output','Filtered')

        subplot(2,1,2)
        hold on
        plot(x,p1,'r')
        plot(x,p2,'b')
        axis([0,250,-pi,pi])

        %Adding highlights to plot
        yl = ylim;
        xl = xlim;

        xBox = [freqRemove(index,1) freqRemove(index,3);
                freqRemove(index,1) freqRemove(index,3);
                freqRemove(index,2) freqRemove(index,4);
                freqRemove(index,2) freqRemove(index,4)];

        yBox = [yl(2) yl(2);
                yl(1) yl(1);
                yl(1) yl(1);
                yl(2) yl(2)];
        patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.1);
        xlabel('\textbf{Frequency [Hz]}')
        ylabel('\textbf{Phase [degrees]}')
        set(gca,'fontname','times')
        set(gca,'Fontsize',12);
        grid on
    end
    

    %% Time domain
    if PlotTime == 1

        figure();
        x = 0:dx:dx*Ntraces-dx;
        t = 0:dt:dt*Nsamples-dt;

        subplot(1,2,1);
        ax = imagesc(x,t,DataY);

        set(0,'defaulttextinterpreter','latex')
        title('\textbf{Objective Gather}')
        xlabel('\textbf{Offset [m]}')
        ylabel('\textbf{Time [s]}')
        colorbar;
        colormap('gray')
        caxis([0 1])
        set(gca,'XAxisLocation','top')
        set(gca,'fontname','times')
        set(gca,'Fontsize',12)

        subplot(1,2,2);
        ax = imagesc(x,t,DataZ);
        title('\textbf{NN Output Gather}')
        xlabel('\textbf{Offset [m]}')
        ylabel('\textbf{Time [s]}')
        colorbar;
        colormap('gray')
        caxis([0 1])
        set(gca,'XAxisLocation','top')
        set(gca,'fontname','times')
        set(gca,'Fontsize',12)

        figure;
        subplot(1,2,1)
        ax = imagesc(x,t,abs(DataY - DataZ));
        title('\textbf{Error}')
        xlabel('\textbf{Offset [m]}')
        ylabel('\textbf{Time [s]}')
        colorbar;
        colormap(flipud(hot))
        caxis([0 0.1])
        set(gca,'XAxisLocation','top')
        set(gca,'fontname','times')
        set(gca,'Fontsize',12)

        set(0,'defaulttextinterpreter','none')
    end

    %Calculate the RMSE
    if CalcSens == 1
    TimeRMSE(index) = mean(mean(sqrt((DataY-DataZ).^2)));
    FreqRMSE(index) = mean(mean(sqrt((f1-f2).^2)));
    PhaseRMSE(index) = mean(mean(sqrt((p1-p2).^2)));
    end

   
end

%Show the RMSE:
if CalcSens == 1 
    TimeRMSE = mean(TimeRMSE);
    FreqRMSE = mean(FreqRMSE);
    PhaseRMSE = mean(PhaseRMSE);

    fprintf('TimeRMSE = %.10f \n', TimeRMSE);
    fprintf('FreqRMSE = %.10f \n', FreqRMSE);
    fprintf('PhaseRMSE = %.10f \n', PhaseRMSE);
end

