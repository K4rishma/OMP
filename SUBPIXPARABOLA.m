function [vector] = SUBPIXPARABOLA(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = find(~((x <= (size(result_conv,2)-1)) & (y <= (size(result_conv,1)-1)) & (x >= 2) & (y >= 2)));
    x(xi) = [];
    y(xi) = [];
    z(xi) = [];
    xmax = size(result_conv, 2);
    vector = NaN(size(result_conv,3), 2, class(result_conv));
    if(numel(x)~=0)
        ip = sub2ind(size(result_conv), y, x, z);
        f0 = result_conv(ip);
        f1 = result_conv(ip-1);
        f2 = result_conv(ip+1);
        peaky = y + (f1-f2)./(2*f1-4*f0+2*f2);
        f0 = result_conv(ip);
        f1 = result_conv(ip-xmax);
        f2 = result_conv(ip+xmax);
        peakx = x + (f1-f2)./(2*f1-4*f0+2*f2);

        SubpixelX=peakx-(interrogationarea(2)/2)-SubPixOffset;
        SubpixelY=peaky-(interrogationarea(1)/2)-SubPixOffset;
        vector(z, :) = [SubpixelX, SubpixelY];  
    end
end

