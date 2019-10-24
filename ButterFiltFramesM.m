function [Frames] = ButterFiltFramesM(Frames, lims, centerFreq, c, FR, ax)
%ButterFiltFrames Filters frames using a zero-phase butterworth filter for
%slow time filtering of data
%this uses a file exchange filter for increased speed
%INPUT:
%   Frames: Data to be filtered - can be any dimension
%   lims: [lower lim, upper lim] velocity to keep in m/s
%   centerFreq: center frequency of pulse (Hz)
%   c: Speed of sound in medium
%   FR: frame rate in Hz
%   ax: axis to perform filtering
%OUTPUT:
%   Frames: Filtered Data
    tic
    low_lim_hz = lims(1) * (2 * centerFreq / c);  % doppler freq Hz
    up_lim_hz = lims(2) * (2 * centerFreq / c);  % doppler freq Hz
    if up_lim_hz >= FR/2
        [b,a] = butter(2, low_lim_hz/(FR/2), 'high');
    else
        [b,a] = butter(2, [low_lim_hz, up_lim_hz]./(FR/2));
    end
    Frames = FiltFiltM(b, a, Frames, ax);
    timer = toc;
    fprintf(['Filtering data took ' num2str(timer,2) 's\n']);

end

