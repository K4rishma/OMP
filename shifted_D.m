function D = shifted_D(im_frame, block_length)
% Creating Dictionary D


%block_length = 16;
kernels = size(im_frame,3);
im_size = block_length.^2;
D = zeros(im_size,im_size,kernels); % Each column, D(:,i) is one of the displacements
shifted_im = flipud(fliplr(im_frame));
% Time complexity is block_length^2 = 256
% Since no of slice in each direction is block_length+1
temp = zeros(block_length , block_length, kernels);

for j = 1:block_length
    for i = 1:block_length
%    for j = 7:10
%     for i = 7:10
        temp(1:i,1:j,:) = shifted_im(1:i,1:j,:);
       
         D(:,i + block_length*(j-1),:) = reshape( flipud(fliplr(temp)),[],kernels);
    end
end

end