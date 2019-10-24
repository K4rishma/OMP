function [vector,maxVals] = correlation2(image1, image2, sw, masked_indices)
% correlation 2 is with shifted D ans see its comparison with
% karishma_Dcorrelation


b_size = size(image1, 1);
kernels = size(image1,4);
n_avg = size(image1,3);
x_disp = zeros(kernels,1);
y_disp = zeros(kernels,1);
max_value = zeros(kernels,1);

meanIm1 = sum(sum(image1))./b_size.^2;
image1 = image1;

D = karishma_Dcorrelation(image2, b_size);
y = reshape(image1, [], 1, n_avg, kernels);
corrmap = zeros(b_size + 1, b_size + 1,kernels);

masked_kernels = 1:kernels;
masked_kernels(masked_indices) = [];

for k = masked_kernels(1) : masked_kernels(end)
cor = zeros(b_size + 1, b_size + 1);

for n = 1 : n_avg
    im2 = permute(D(:,:,n,k), [2,1,3,4]); 
    %temp = ((D(:,:,n,k)./std(D(:,:,n,k))).')*y(:,1,n,k)./std2(image1(:,:,n,k));    
    %ver = xcorr2(image2(:,:,n,k)./std2(image2(:,:,n,k)),image1(:,:,n,k)./std2(image1(:,:,n,k)));
    Im1 = ( y(:,1,n,k) - mean(y(:,1,n,k)) )./std(y(:,1,n,k));
    temp = ((D(:,:,n,k)).')*Im1;  
    temp = reshape(temp, b_size+1, b_size+1);  
    temp = temp./(size(image1,1)*size(image1,2));
    cor = cor + temp;
end


cor = (cor./n_avg);
cor = sw(:,:,k).*cor;
corrmap(:,:,k) = cor;
% max_value(k) = max(max(cor));
% [y_disp(k), x_disp(k)] = find(cor == max_value(k) )   ; % image 1 and 2 changed?, see how correlation works (abs?)
% 
 end
% 
% vert = y_disp;
% hort = x_disp;
% 
% 
% 
 ksize = Opts.kSize;
 corrmap(:,:,masked_indices) = NaN;
 [vector, maxVals] = subpix_fit(corrmap, ksize);
% 
% 
% Y_disp = -(b_size)/2 - 1 + vert ; %see this || 
% X_disp = -(b_size)/2  - 1 + hort ;
% 
% vector = zeros(kernels,2);
% vector = [(X_disp(:)), (Y_disp(:))];


end