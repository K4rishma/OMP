function [polyline] = drawPoly(im,drange)
    if (nargin==1)
        drange=[min(im(:)),max(im(:))];
    end
    scnsize = get(0,'ScreenSize');pos=[10 scnsize(4)/2 scnsize(3)/2 scnsize(4)/2]; %window position on screen
    disp('Draw the polyline you want to measure. Do not close the polygon');
    hfig1=figure('position', pos); hAx = axes; imshow(im,drange,'InitialMagnification',100,'Parent', hAx);
    %pick rectangle region in source frame for tracking:
    polyf = impoly;
    polyline = polyf.getPosition;
    close(hfig1); 

end