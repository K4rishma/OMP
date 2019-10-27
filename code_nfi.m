close all;
clear all;

%% Parameter Definitions
% Objective: To detect rising and falling edges and smoothen the measured signal
% % measured signal
% meas_signal          % measured signal
% f1 f2 f3             % Hz, various frequencies present in the measured signal
% A1 A2 A3             % amplitude of above frequencies
% fr                   % Hz, reference frequency of the measured signal
% time_stamps          % sec, time observations for the measured signal

% % recovered signals
% recov_sig            % recovered signal : filtered measured signal using moving average
% pre                  % number of previous samples used in moving average
% variation            % standard deviation of moving average
% lse                  % least square error between recovered and measured signal
% threshold            % value for peak detection in lse
% points               % Array containing the position of mid point of rising and falling edges 
% buffer               % number of indices near the points that will be used for smoothening purpose 
% recov_sig2           % another recovered signal that does not smoothen the edges
%% Improvements
% measured signal can have varying refrence frequencies instead of only one
% recov_sig2 can be smoothly tuned, near to the rising and falling edges (see the visualisation) 
% appication-wise, threshold and buffer can be made adaptive to the system 
%% Input Parameters
f1 = 50; f2 = 120; f3 = 200;
A1 = 0.5; A2 = 0.5; A3 = 0.5;
time_stamps = (1: 0.2: 200)./200;
fr = 10;

pre = 10;
threshold = 18;
buffer = 2;
%% Measured signal simulation
N = length(time_stamps);
sin_r = sin(2*pi*fr*time_stamps);
ref_signal = 10.*(sin_r > 0) + 5.*(sin_r <= 0);
fluctuations = A1.*sin(2*pi*f1*time_stamps) + A2.*sin(2*pi*f2*time_stamps) + A3.*sin(2*pi*f3*time_stamps);
meas_signal = ref_signal + fluctuations + randn(1,N);

%% Signal recovery
variation = zeros(1,N);

recov_sig = meas_signal;
recov_sig2 = meas_signal;

for  i = 1 : N - pre
    recov_sig(pre + i - 1 ) = mean( meas_signal(i : pre + i - 1 ));
    variation(pre + i - 1 ) = std( meas_signal(i : pre + i - 1 ));
end

lse = (meas_signal - recov_sig).^2 + variation.^2;

detect = double((lse > threshold));
differece = detect;
differece(2:end) = detect(2:end) - detect(1:end-1);
points = (differece == 1);
index = find(points);                                                       % storing indices at which lse just crosses the threshold


for i = 2:length(index)
    
    if (index(i) - index(i-1)) >= (pre + 10)
        signal_matrix = [];
        start = index(i-1) + buffer;
        finish = index(i) - buffer;
        inc = finish - start + 1 ;
        for j = pre:inc
            column = j - pre + 1;
            signal_matrix(:,column) = meas_signal(start + column - 1 : ...
                start + column - 1 + pre - 1).';
        end
        recov_sig2( start + pre - 1 : finish) = mean(signal_matrix);
    else
        index(i) = index(i-1);
    end
    
end

points_m = zeros(1,N);
points_m(index) = 1;                                                        % moidified points containing one point per edge

%% Visualisation


figure;
subplot(2,1,1);
plot(ref_signal); 
hold on 
plot(7.5 + sin_r);
hold off
legend('square','sine')
subplot(2,1,2);
plot(meas_signal);
legend('measured')
xlabel('time')


figure;
plot(meas_signal);
hold on
plot(recov_sig);
plot(recov_sig2);
plot(points_m);
hold off;
legend('measured','recov','recov2','edge detection');
xlabel('time'); ylabel('Intensity');

figure;
plot(meas_signal);
hold on
plot(recov_sig2);
hold off;
legend('measured','recov2');
xlabel('time'); ylabel('Intensity');