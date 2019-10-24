function [idx] = live_view_quiv(image, x, y, u, v, sf, idx, drange)
%live_view_quiv Plots quiver for a frame over the mean of the image stack
%used to calculate the quiver
%INPUT
%image: NxMxEnsembles image frames (envelope)
%x:     nxm x coordinates of vector bases
%y:     nxm y coordinates of vector bases
%u:     nxm x displacements of vectors
%v:     nxm y displacements of vectors
%sf:    scale factor for displacement vectors (larger = larger vectors)
%idx:   a counter for counting through the frames. Should init with 0
%drange: dynamic range to show images at eg. [-50 0]
%OUTPUT
%idx:   return the counter so it can be used in the next iteration
persistent im quiv norm
if idx==0
       figure();
       meanim = mean(image,3);
       norm = max(meanim(:));
       im = imagesc(20*log10(meanim./norm), drange);
       colormap gray
       hold on
       quiv = quiver(x,y,u.*sf,v.*sf,0);
       idx=1;
       title(idx)
       drawnow;
else
    meanim = mean(image,3);
    norm = max([norm, max(meanim(:))]);
    im.CData = 20*log10(meanim./norm);
    quiv.UData = u.*sf;
    quiv.VData = v.*sf;
    idx = idx+1;
    title(idx)
    drawnow;
end

end

