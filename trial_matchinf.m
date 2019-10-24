clc;
close all; clear all;

% NCC 
b1 = randn(16); b2 = randn(16);
b1(7:10,7:10) = 10; b2(1:4,1:4) = 10;

m1 = mean(mean(b1));   m2 = mean(mean(b2));
forier = ifft2((conj(fft2( b1 - m1)).*(fft2( b2 - m2))));
num = fftshift(fftshift(real(forier), 1), 2);
den = std2(b1)*std2(b2)*256; % size missing

ncc = num./den;
[r1, c1 ] = find( ncc == max(max(ncc)), 1 );
r1 = r1-1 - 8
c1 = c1-1 - 8
% correlation
im1 = (b1 - m1)./std2(b1);
im2 = (b2 - m2)./std2(b2);

corr = xcorr2(im2, im1);
%corr = corr./256;
[r2, c2 ] = find( corr == max(max(corr)), 1 );
r2 = r2 - 16
c2 = c2 - 16