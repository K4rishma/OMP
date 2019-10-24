function [utable, vtable, removed] = median_outlier(utable,vtable,epsilon,thresh,b)
% from Westerweel, J. and Scarano, F., “Universal outlier detection for PIV data,” 
% Exp. Fluids, vol. 39, no. 6, pp. 1096–1100, Dec. 2005.
    [J,I]=size(utable);
    normfluct=zeros(J,I,2, class(utable));
    for c=1:2
        if c==1
            velcomp=utable;
        else
            velcomp=vtable;
        end
        clear neigh
        for ii = -b:b
            for jj = -b:b
                neigh(:, :, ii+2*b, jj+2*b)=velcomp((1+b:end-b)+ii, (1+b:end-b)+jj);
            end
        end
        neighcol = reshape(neigh, size(neigh,1), size(neigh,2), (2*b+1)^2);
        neighcol2= neighcol(:,:, [(1:(2*b+1)*b+b) ((2*b+1)*b+b+2:(2*b+1)^2)]);
        neighcol2 = permute(neighcol2, [3, 1, 2]);
        med=median(neighcol2);
        velcomp = velcomp((1+b:end-b), (1+b:end-b));
        fluct=velcomp-permute(med, [2 3 1]);
        res=neighcol2-repmat(med, [(2*b+1)^2-1, 1,1]);
        medianres=permute(median(abs(res)), [2 3 1]);
        normfluct((1+b:end-b), (1+b:end-b), c)=abs(fluct./(medianres+epsilon));
    end
    info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
    removed = sum(info1(:));
    utable(info1==1)=NaN;
    vtable(info1==1)=NaN;
end