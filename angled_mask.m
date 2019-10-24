function mask = angled_mask(hor_centre,vert_centre,theta_range,kgrid,radius,piv)
xCenter = hor_centre;
yCenter = vert_centre;
theta = theta_range(1) : 4 :theta_range(2);
radius = radius + 8;
x = radius * cosd(theta) + xCenter;
y = radius * sind(theta) + yCenter;
x = [x,xCenter];
y = [y,yCenter];
% plot(x, y);
% axis square;
% xlim([0 20]);
% ylim([0 20]);
% grid on;
mask_temp = poly2mask(x,y,kgrid.Ny, kgrid.Nx);
ker_dim1 = size(piv.U,1);
ker_dim2 = size(piv.U,2);
% mask = imresize(mask,[ker_dim1,ker_dim2]);
ss1 = piv.ss1s{1};
ss1_mask = ss1(:,:,1,:);
seg_mask = mask_temp(ss1_mask);
kernels = size(ss1_mask,4);
seg_mask = reshape(seg_mask,[],kernels);
mean_mask = nanmean(seg_mask);
mask = reshape(mean_mask,[ker_dim1,ker_dim2]);
mask(mask>0) = 1;
mask(mask==0)= NaN;
% Mask the image using bsxfun() function
% maskedRgbImage = bsxfun(@times, rgbImage, cast(mask, class(rgbImage)));