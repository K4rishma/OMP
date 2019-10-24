function [Frames] = ButterFiltFrames( Frames, lims, centerFreq, c, FR , order)
%ButterFiltFrames Filters frames using a zero-phase butterworth filter for
%slow time filtering of data
%INPUT:
%   Frames: Data to be filtered - can be any dimension
%   lims: [lower lim, upper lim] velocity to keep in m/s
%   centerFreq: center frequency of pulse (Hz)
%   c: Speed of sound in medium
%   FR: frame rate in Hz
%   order: filter order
%OUTPUT:
%   Frames: Filtered Data
    tic
    Frames = permute(Frames, [4,1,2,3]);  % frames, r, theta, angles
    low_lim_hz = 2 * centerFreq * lims(1) / c;  % doppler freq Hz
    up_lim_hz = 2 * centerFreq * lims(2) / c;  % doppler freq Hz
    fprintf('Cutoff low = %4.2f Hz, cutoff high = %4.2f Hz, Nyquist = %4.2f Hz', low_lim_hz, up_lim_hz, FR/2);
    if up_lim_hz >= FR/2
        [z,p,k] = butter(order, low_lim_hz/(FR/2), 'high');
    else
        [z,p,k] = butter(order, [low_lim_hz, up_lim_hz]./(FR/2));
    end
    [sos, g] = zp2sos(z,p,k);
    Frames = single(filtfilt(sos,g,double(Frames)));
    Frames = permute(Frames, [2,3,4,1]); % r, theta, angles, frames
    timer = toc;
    fprintf(['Filtering data took ' num2str(timer,2) 's\n']);

end

