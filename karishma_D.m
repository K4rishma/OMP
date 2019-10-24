function D = karishma_D(im_frame, block_length)
% Creating Dictionary D


%block_length = 16;
kernels = size(im_frame,3);
im_size = block_length.^2;
D = zeros(im_size,block_length.^2,kernels); % Each column, D(:,i) is one of the displacements

% Time complexity is block_length^2 = 256
% Since no of slice in each direction is block_length+1
for j = 1:block_length
    for i = 1:block_length
        D(:,i+block_length*(j-1),:) = reshape(im_frame(i:i+block_length-1,...
            j:j+block_length-1,:),[],kernels);
    end
end

end



