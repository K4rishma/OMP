%% Plot spectrum
numfft= 128;
winsize = 32;
ovrlp = 8;
inmask = find(imMask);
framesInMask = reshape(Frames, [], P.na, size(Frames,3)/P.na);
framesInMask = framesInMask(inmask,:,:);
dims = size(spectrogram(squeeze(framesInMask(1,1,:)), winsize, ovrlp, numfft, FR));
fspec = zeros([dims, size(framesInMask,1), P.na]);
for j = 1:P.na
    for i = 1:size(framesInMask,1)
        fspec(:,:,i,j) = spectrogram(squeeze(framesInMask(i,j,:)), winsize, ovrlp, numfft, FR);
    end
end
temp = zeros([size(fspec,1), size(fspec,2),size(Frames,1)*size(Frames,2), P.na]);
temp(:,:,inmask,:) = fspec;
temp = reshape(temp, [size(fspec,1), size(fspec,2),size(Frames,1),size(Frames,2), P.na]);    