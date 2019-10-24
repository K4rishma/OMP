 function [corr] = normxcorr(image1, image2)
    meanIm1 = mean(mean(mean(image1)));
    meanIm2 = mean(mean(mean(image2)));
    
    padding = zeros(size(image2));
    padding(1:size(image1,1),1:size(image1,1),:,:) = image1 - meanIm1; 
    
    stdIm1 = permute(std(reshape(image1,[], size(image1,4))),[3,4,1,2]);
    stdIm2 = permute(std(reshape(image2,[], size(image2,4))),[3,4,1,2]);
    %corrFFT = conj(fft2(image1-meanIm1, size(image2,1), size(image2,2))).*fft2(image2-meanIm2);
    corrFFT = conj(fft2(padding)).*fft2(image2-meanIm2);
    corrFFT = mean(corrFFT,3);  % mean can be done here saving time (check for nans!!!)
    result_conv = fftshift(fftshift(real(ifft2(corrFFT)), 1), 2);
    normFact = (size(image1,1)*size(image1,2)).*stdIm1.*stdIm2;
%     corr = result_conv(8:23,8:23,:,:)./normFact;
    corr = result_conv(1:size(image1,1),1:size(image1,1),:,:)./normFact;


% meanIm1 = mean(mean(mean(image1)));
%     meanIm2 = mean(mean(mean(image2)));
%     stdIm1 = permute(std(reshape(image1,[], size(image1,4))),[3,4,1,2]);
%     stdIm2 = permute(std(reshape(image2,[], size(image2,4))),[3,4,1,2]);
%     corrFFT = conj(fft2(image1-meanIm1)).*fft2(image2-meanIm2);
%     corrFFT = mean(corrFFT,3);  % mean can be done here saving time (check for nans!!!)
%     result_conv = fftshift(fftshift(real(ifft2(corrFFT)), 1), 2);
%     normFact = (size(image1,1)*size(image1,2)).*stdIm1.*stdIm2;
%     corr = result_conv./normFact;
end