function [indxs] = kernalize(mini, maxi,step, nK, dDims, ksize)
    % creates tensor of linear indices for referencing frame data
    % to allow for vectorized operation (no matlab loops when
    % processing xcorr)
    s0y = permute(repmat((mini(1):step(1):maxi(1))-1,[nK(2),1,nK(3)]),[2 1 3]);
    s0x = permute(repmat(((mini(2):step(2):maxi(2))-1).*dDims(1),[nK(1),1,nK(3)]),[1 2 3]);
    s0t = permute(repmat(((mini(3):step(3):maxi(3))-1).*(dDims(1)*dDims(2)),[nK(1),1,nK(2)]),[1 3 2]);
    s0 = s0x+s0y+s0t;
    s0 = uint32(reshape(s0,[],ksize(3))); %all indexes for top left corners of kernels
    s1y = permute(repmat(1:ksize(1),[ksize(2),1,ksize(3)]),[2,1,3]);
    s1x = permute(repmat(((1:ksize(2))-1).*dDims(1),[ksize(1),1,ksize(3)]),[1,2,3]);
    %     s1t = permute(repmat(((1:ksize(3))-1).*(dDims(1)*dDims(2)),[ksize(1),1,ksize(2)]),[1,3,2]);
    s1 = uint32(s1y+s1x); % all indexes inside a kernel
    indxs = repmat(s1,[1, 1, 1, size(s0,1)])+repmat(permute(s0,[4,3,2,1]),[ksize(1),ksize(2),1,1]); %indexes of each pixel in each kernel
end

