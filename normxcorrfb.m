function [corr] = normxcorrfb(image1, image2)
    % assumes even kernel size. so circular shift is required by 1 unit in
    % x/y. For odd kernel size no circular shift is necessary and should be
    % added as an input.
    meanIm1 = mean(mean(image1));
    meanIm2 = mean(mean(image2));
    stdIm1 = permute(std(reshape(image1,[],size(image1,3), size(image1,4))),[4,1,2,3]);
    stdIm2 = permute(std(reshape(image2,[],size(image2,3), size(image2,4))),[4,1,2,3]);
    corrFFT = conj(fft2(image1-meanIm1)).*fft2(image2-meanIm2);
    result_convfw = fftshift(fftshift(real(ifft2(corrFFT)),1),2);
    corrFFT = conj(fft2(image2-meanIm2)).*fft2(image1-meanIm1);
    result_convbw = fftshift(fftshift(real(ifft2(corrFFT)),1),2);
    result_conv = (result_convfw + circshift(flip(flip(result_convbw,1),2),[1,1,0,0]))./2;
    normFact = (size(image1,1)*size(image1,2)).*stdIm1.*stdIm2;
    corr = result_conv./normFact;
    
    %debug
%     subplot(1,3,1)
%     imagesc(result_convfw(:,:,1,1))
%     subplot(1,3,2)
%     imagesc(circshift(flip(flip(result_convbw(:,:,1,1),1),2),[1,1,0,0]))
%     subplot(1,3,3)
%     imagesc(result_conv(:,:,1,1))
    
    
end