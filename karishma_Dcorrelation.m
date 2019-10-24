function D = karishma_Dcorrelation(im_frame, block_length)
% Creating Dictionary D
% difference from karishma_D is image dimensions are different 

%block_length = 16;
kernels = size(im_frame,4);
n = size(im_frame,3);
im_size = block_length.^2;
D = zeros(im_size, (block_length +1).^2, n, kernels); % Each column, D(:,i) is one of the displacements

% Time complexity is block_length^2 = 256
% Since no of slice in each direction is block_length+1
for j = 1:block_length + 1
    for i = 1:block_length + 1
        clear L
        L = im_frame(i:i+block_length-1,j:j+block_length-1,:,:);
%         L = L - mean( sum(sum(L))./im_size );
%         D(:,i+(block_length +1)*(j-1),:,:) = reshape(L,[],n,kernels);
          L = reshape(L,[],n,kernels);
          L = (L - mean(L))./std(L);
          D(:,i+(block_length +1)*(j-1),:,:) = reshape(L,[],n,kernels);
    end
end

end