function [Frames] = FIRFiltFrames( Frames, lims, centerFreq, c, FR, order )
%FIRFiltFrames Filters frames using a linear-phase Chebyshev filter for
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
    low_lim_hz = lims(1) * (2 * centerFreq / c);  % doppler freq Hz
    up_lim_hz = lims(2) * (2 * centerFreq / c);  % doppler freq Hz
    if up_lim_hz >= FR/2
        b = fir1(order,low_lim_hz/(FR/2),'high',chebwin(order+1,45));
    else
        b = fir1(order, [low_lim_hz, up_lim_hz]./(FR/2),chebwin(Opts.filterOrder+1,45));
    end
    Frames = filter(b,1,Frames);
    Frames = permute(Frames, [2,3,4,1]); % r, theta, angles, frames
    timer = toc;
    fprintf(['Filtering data took ' num2str(timer,2) 's\n']);

end

