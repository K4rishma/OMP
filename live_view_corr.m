function [idx] = live_view_corr(corrim, u, v, sf, idx)
%live_view_quiv Plots quiver for a frame over the mean of the image stack
%used to calculate the quiver
%INPUT
%corrim:nxm correlation map image frames (envelope)
%u:     nxm x displacements of vectors
%v:     nxm y displacements of vectors
%sf:    scale factor for displacement vectors (larger = larger vectors)
%idx:   a counter for counting through the frames. Should init with 0
%OUTPUT
%idx:   return the counter so it can be used in the next iteration
persistent im quiv
if idx==0
       figure();
       im = imagesc(corrim, [0,1]);
       colormap hot
       hold on
       quiv = quiver(u.*sf,v.*sf,0);
       idx=1;
       title(idx)
       drawnow;
else

    im.CData = corrim;
    quiv.UData = u.*sf;
    quiv.VData = v.*sf;
    idx = idx+1;
    title(idx)
    drawnow;
end

end

