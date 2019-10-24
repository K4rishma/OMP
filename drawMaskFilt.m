function [imMask, ROI, poly] = drawMaskFilt(dataIQ, hp, tax, ang_ax)
    if nargin<3
        tax = 4; % time dimension
        ang_ax = 3;  % dimension with angles
    end
    if nargin<2
        hp = 0.5;
    end
    [b, a] = butter(4, hp, 'high');
    im = filter(b,a,dataIQ,[], tax);  % dont care about 0-phase, just a quick filter for seeing image.  
    im = sum(abs(im),ang_ax); % incoherent compounding
    im = max(im, [], tax);
    im = squeeze(20*log10(im));
    [imMask, ROI, poly] = drawMask(im);

end