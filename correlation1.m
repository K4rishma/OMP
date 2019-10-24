function [vector,hort,vert] = correlation1(image1, image2, sw, masked_indices)
%image1 and image2 same size xcorr shows result

meanIm1 = mean(mean(mean(image1)));
meanIm2 = mean(mean(mean(image2)));

b_size = size(image1, 1);
se_size = size(image2, 1);
kernels = size(image1,4);
n_avg = size(image1,3);
x_disp = zeros(kernels,1);
y_disp = zeros(kernels,1);

image1 = image1 - meanIm1;
image2 = image2 - meanIm2;

%D = karishma_Dcorrelation(image2, b_size)

for k = 1:kernels
    cor = zeros(2*b_size-1);
    for n = 1:n_avg
       % image2_D = D(:,)
        temp = xcorr2(image2(:,:,n,k)./std2(image2(:,:,n,k)),image1(:,:,n,k)./std2(image1(:,:,n,k)));
        
        cor = cor + temp;
    end
    cor = cor./n_avg;
    [y_disp(k),x_disp(k)] = find(cor == max(max(cor)) )   ;
end
vert = y_disp;
hort = x_disp;

Y_disp = -b_size + vert ;
X_disp = -b_size + hort ;

vector = zeros(kernels,2);
vector = [X_disp(:),Y_disp(:)];


end