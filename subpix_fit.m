function [vector, maxVals] = subpix_fit(corrmap, ksize)
%subpix_fit Find the subpixel displacement by fitting a function over the
%correlation map. Function is chosen in the Opts class.
%INPUT
%   corrmap - 2D set of correlation maps. (MxNxKernels)
%   ksize - kernel size for this iteration
%OUTPUT
%   vector - 2xKernels array of displacements (y,x)
%   maxVals - maximum xcorr values for use later
    if (rem(ksize,2) == 0) %for the subpixel displacement measurement
        SubPixOffset=1;
    else
        SubPixOffset=0.5;
    end
    [maxVals,maxInds] = nanmax(reshape(corrmap,[],size(corrmap,3)),1);
    maxInds = maxInds';
   
%  kernels = size(corrmap,3);   
%  b_size = size(corrmap,1);
% [vert,hort] = ind2sub([b_size,b_size],maxInds);
% Y_disp = vert - b_size/2 - 1;
% X_disp = hort - b_size/2 -1;
% 
% 
% vector = zeros(kernels,2);
% vector = [X_disp(:),Y_disp(:)];
% vector(masked_indices,:) = NaN;
    
    
    [y, x, ~] = ind2sub(size(corrmap), maxInds);
    z = [1:size(corrmap,3)]';
    [z1, zi] = sort(z);
    % we need only one peak from each couple pictures
    dz1 = [z1(1); diff(z1)];
    i0 = find(dz1~=0);
    x1 = x(zi(i0));
    y1 = y(zi(i0));
    z1 = z(zi(i0));
   
    if Opts.subPixFinder==1
        [vector] = SUBPIXGAUSS (corrmap,ksize, x1, y1, z1,SubPixOffset);
    elseif Opts.subPixFinder==2
        [vector] = SUBPIX2DGAUSS (corrmap,ksize, x1, y1, z1,SubPixOffset);
    elseif Opts.subPixFinder==3
        [vector] = SUBPIXCENTROID (corrmap, ksize, x1, y1, z1, SubPixOffset);
    elseif Opts.subPixFinder==4
        [vector] = SUBPIXPARABOLA (corrmap, ksize, x1, y1, z1, SubPixOffset);
    elseif Opts.subPixFinder==5
        [vector] = SUBPIX2DPARAB (corrmap, ksize, x1, y1, z1, SubPixOffset);
    end

end

