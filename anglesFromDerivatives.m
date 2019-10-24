function angles = anglesFromDerivatives(acquisitions,mask)
%figure; imagesc(mean(abs(BF_filtered),3));title('filtered PDI');
I = (mean(abs(acquisitions(:,:,1,:)),4));
%d_mask = I > graythresh(I);%- 0.60*graythresh(I);
d_mask = mask;
masked_I = I.*d_mask;
%figure; imagesc(masked_I);

dx_masked = 0.*d_mask;
dy_masked = 0.*d_mask;

dx_masked(:,2:end) = diff(masked_I,1,2);
dy_masked(2:end,:) = diff(masked_I,1,1);
dx_scaled = dx_masked(:,:);
dy_scaled = dy_masked(:,:);
%figure; imagesc([dx_scaled dy_scaled]);
%Opts.dsize= 1;
smoothner = ones(Opts.dsize) ;
dx_smooth = conv2(dx_scaled,smoothner,'same'); dy_smooth = conv2(dy_scaled,smoothner,'same');

%figure; imagesc([dx_smooth dy_smooth]);

angles = atand((dx_smooth)./(dy_smooth));
angles = medfilt2(angles,[Opts.dsize Opts.dsize]);

%find(angles(i,:),1)

% Rigged
angles = -90.*d_mask;
angles(angles==0) = NaN;

%  figure; imagesc(angles);
%  figure; imagesc(d_mask);
%  figure; imagesc(I);
end