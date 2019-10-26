function SR_indices = Indices_Dcorrelation(sr_size, block_length)
% Creating indices for Dictionary D
% sr_size is the cut size for image 2 [block_length*2 block_length*2 n kernels ]
% created to optimise timing

kernels = sr_size(4);
n = sr_size(3);
im_vector = block_length.^2;
SR_indices = zeros(im_vector, (block_length +1).^2, n, kernels);

im_frame = zeros(sr_size);

for j = 1:block_length + 1
    for i = 1:block_length + 1
        im_frame = 0.*im_frame;
        clear indices
        im_frame(i:i+block_length-1,j:j+block_length-1,:,:) = 1;
        indices = find(im_frame);
        indices = reshape(indices,[],n,kernels);
        SR_indices(:,i+(block_length +1)*(j-1),:,:) = indices;
    end
end

end