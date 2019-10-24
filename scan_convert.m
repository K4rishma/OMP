function [FramesCart] = scan_convert(FramesPolar, R, THETA, Z, X)
%scan_convert - scan converts FramesPolar on a grid defined by res and
%origin to grid.
%   FramesPolar: Polar domain beamformed data
%   R: ndgrid of radial positions of each pixel in FramesPolar
%   THETA = ndgrid of theta positions of each pixel in FramesPolar
%   Z = ndgrid of depth positions to scan convert - same units as R
%   X = ndgrid of lateral positions to scan convert - same units as R
%% Scan convert the fast way
    R_ = sqrt(X.^2 + Z.^2);
    THETA_ = atan2(X,Z);

    fint = griddedInterpolant(R,THETA, abs(FramesPolar(:,:,1)),'linear');
    fint.ExtrapolationMethod = 'none';
    FramesCart = zeros([size(Z) size(FramesPolar,3)],'single');
    for i = 1:size(FramesPolar,3)
        fint.Values = FramesPolar(:,:,i);
        FramesCart(:,:,i) = fint(R_,THETA_);
    end
end