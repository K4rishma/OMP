function [vector,hort,vert] = correlation1(image1, image2, sw, masked_indices) 
    meanIm1 = mean(mean(mean(image1)));
    meanIm2 = mean(mean(mean(image2)));

    b_size = size(image1, 1);
    se_size = size(image2, 1);
    kernels = size(image1,4);
    n_avg = size(image1,3);
x_disp = zeros(kernels,1); y_disp = zeros(kernels,1);
for k = 1:kernels
    cor = [];
    for n = 1:n_avg
temp = xcorr2(image2(:,:,k)./std2(image2(:,:,k)),image1(:,:,k)./std2(image1(:,:,k)));
% cor = xcorr2(image1(:,:,k),image2(:,:,k));
cor = cor + temp;
    end
    corr = corr./n_avg;
[y_disp(k),x_disp(k)] = find(cor == max(max(cor)) )   ;
end
vert = y_disp;
hort = x_disp;
% Y_disp = nanmean(-b_size + vert - 1 ,2);
% X_disp = nanmean((-b_size + hort - 1 ),2);
Y_disp = -b_size + vert ;
X_disp = -b_size + hort ;

vector = zeros(kernels,2);
vector = [X_disp(:),Y_disp(:)];
%vector(masked_indices,:) = NaN;

end