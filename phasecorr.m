function [corr] = phasecorr(image1, image2)
    corrFFT = conj(fft2(image1)).*fft2(image2);
    corrFFT = corrFFT./abs(corrFFT);
    corrFFT = mean(corrFFT,3);  % mean can be done here saving time (check for nans!!!)
    corr = fftshift(fftshift(real(ifft2(corrFFT)), 1), 2);
end