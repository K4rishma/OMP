function [vector] = SUBPIX2DGAUSS(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = find(~((x <= (size(result_conv,2)-1)) & (y <= (size(result_conv,1)-1)) & (x >= 2) & (y >= 2)));
    x(xi) = [];
    y(xi) = [];
    z(xi) = [];
    xmax = size(result_conv, 2);
    vector = NaN(size(result_conv,3), 2, class(result_conv));
    if(numel(x)~=0)
        c10 = zeros(3,3, length(z), class(result_conv));
        c01 = c10;
        c11 = c10;
        c20 = c10;
        c02 = c10;
        ip = sub2ind(size(result_conv), y, x, z);

        for i = -1:1
            for j = -1:1
                %following 15 lines based on
                %H. Nobach � M. Honkanen (2005)
                %Two-dimensional Gaussian regression for sub-pixel displacement
                %estimation in particle image velocimetry or particle position
                %estimation in particle tracking velocimetry
                %Experiments in Fluids (2005) 38: 511�515
                c10(j+2,i+2, :) = i*log(result_conv(ip+xmax*i+j));
                c01(j+2,i+2, :) = j*log(result_conv(ip+xmax*i+j));
                c11(j+2,i+2, :) = i*j*log(result_conv(ip+xmax*i+j));
                c20(j+2,i+2, :) = (3*i^2-2)*log(result_conv(ip+xmax*i+j));
                c02(j+2,i+2, :) = (3*j^2-2)*log(result_conv(ip+xmax*i+j));
                %c00(j+2,i+2)=(5-3*i^2-3*j^2)*log(result_conv_norm(maxY+j, maxX+i));
            end
        end
        c10 = (1/6)*sum(sum(c10));
        c01 = (1/6)*sum(sum(c01));
        c11 = (1/4)*sum(sum(c11));
        c20 = (1/6)*sum(sum(c20));
        c02 = (1/6)*sum(sum(c02));
        %c00=(1/9)*sum(sum(c00));

        deltax = squeeze((c11.*c01-2*c10.*c02)./(4*c20.*c02-c11.^2));
        deltay = squeeze((c11.*c10-2*c01.*c20)./(4*c20.*c02-c11.^2));
        peakx = x+deltax;
        peaky = y+deltay;

        SubpixelX = peakx-(interrogationarea(2)/2)-SubPixOffset;
        SubpixelY = peaky-(interrogationarea(1)/2)-SubPixOffset;

        vector(z, :) = [SubpixelX, SubpixelY];
    end
end
