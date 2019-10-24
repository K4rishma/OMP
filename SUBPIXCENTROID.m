function [vector] = SUBPIXCENTROID(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = find(~((x <= (size(result_conv,2)-1)) & (y <= (size(result_conv,1)-1)) & (x >= 2) & (y >= 2)));
    x(xi) = [];
    y(xi) = [];
    z(xi) = [];
    xmax = size(result_conv, 2);
    vector = NaN(size(result_conv,3), 2, class(result_conv));
    if(numel(x)~=0)
        ip = sub2ind(size(result_conv), y, x, z);
        fn = (y-1).*result_conv(ip-1)+y.*result_conv(ip)+(y+1).*result_conv(ip+1);
        fd = result_conv(ip-1)+result_conv(ip)+result_conv(ip+1);
        peaky = fn./fd;
        fn = (x-1).*result_conv(ip-xmax)+x.*result_conv(ip)+(x+1).*result_conv(ip+xmax);
        fd = result_conv(ip-xmax)+result_conv(ip)+result_conv(ip+xmax);
        peakx = fn./fd;

        SubpixelX=peakx-(interrogationarea(2)/2)-SubPixOffset;
        SubpixelY=peaky-(interrogationarea(1)/2)-SubPixOffset;
        vector(z, :) = [SubpixelX, SubpixelY];  
    end
end

