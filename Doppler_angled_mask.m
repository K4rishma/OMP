function dmask = Doppler_angled_mask(hor_centre,vert_centre,theta_range,kgrid,radius,piv)
xCenter = hor_centre;
yCenter = vert_centre;
theta = theta_range(1) : 4 :theta_range(2);
radius = radius + 5;
x = radius * cosd(theta) + xCenter;
y = radius * sind(theta) + yCenter;
x = [x,xCenter];
y = [y,yCenter];
% plot(x, y);
% axis square;
% xlim([0 20]);
% ylim([0 20]);
% grid on;
dmask = poly2mask(x,y,kgrid.Ny, kgrid.Nx);
dmask = double(dmask);
dmask(dmask == 0)= NaN;
% Mask the image using bsxfun() function
% maskedRgbImage = bsxfun(@times, rgbImage