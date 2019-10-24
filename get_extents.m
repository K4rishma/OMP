function [mini, maxi, n_el] = get_extents(ksize, step, dDims,cas)
    % finds pixel positions of first and last kernel and number of kernels that fit in
    % ROI
    if cas == 2
    mini=single([1,1,1]); %inner index of first kernel
    n_el = floor((dDims - ksize)./step + 1);
    else
        mini = single([ksize(1)/2+1,ksize(2)/2+1,1]);
        dDims(1:2) = dDims(1:2) - (2*(mini(1:2)-1)); 
        n_el = floor((dDims - ksize)./step + 1);
    end
     
     %number of kernels that fit into the data dimensions
    maxi = (n_el-1).*step+1+(mini - 1); % inner index of last kernel
%     centershift = floor((dDims-(maxi+ksize-1) - mini)./2); %determine amount to move kernels by to shift them to center of data
%     centershift(centershift<0)=0;% negative shift values not allowed (outside data space)
%     mini = mini+centershift;% shift first kernel inner corner
%     maxi = maxi+centershift;% shift last kernel inner corner
end

