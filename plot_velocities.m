 close all hidden
%  clear all
% clc
%% Load the files %%%
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\thinDiscRBC0_Ne_100_4e5_dx_.mat')
% addpath('statsLib');  % this assumes 2DEchoPIV folder is working dir.
% addpath('Plotting');  % this assumes 2DEchoPIV folder is working dir.
% addpath('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\k-wave-toolbox-version-1.2.1\k-Wave');
%% Assign the values
% run PIVforJason.m % 
%% Visualization %%%%
disc_radius = 60;
range = 10;
SimTime = kgrid.dt*kgrid.Nt;
velocity = 3.99;%0.005*20; % maximum blood velocity % to check freq resoulution max blood velocity*num of angles* frame skip(optional)
theta = atand(velocity*SimTime/(kgrid.dx*disc_radius));

V1 = makeDisc(kgrid.Ny, kgrid.Nx, 120, 105, disc_radius);
V2 = makeDisc(kgrid.Ny, kgrid.Nx, 120, 105, disc_radius - range);
V1 = V1 - V2;
actual_velocities = zeros(kgrid.Ny, kgrid.Nx);
for i = 1:range
    V1_temp = makeDisc(kgrid.Ny, kgrid.Nx, 120, 105, disc_radius-(i-1));
    V2_temp = makeDisc(kgrid.Ny, kgrid.Nx, 120, 105, disc_radius-1-(i-1));
    temp_velocity = velocity*(disc_radius-(i-1))/disc_radius;
    V_temp = temp_velocity.*(V1_temp - V2_temp);
    actual_velocities = actual_velocities + V_temp;
end
actual_velocities(actual_velocities == 0) = NaN;

% after post process conversion of displacement to velocities
horz_vel = pivfilt.U;
vert_vel = pivfilt.V;

measured_velocities = sqrt(horz_vel.^2 + vert_vel.^2);
realizations = size(measured_velocities,3);

figure;
max_vel = nanmax(measured_velocities(:));
for i = 1 : realizations
    imagesc(measured_velocities(:,:,i),[0 max_vel]);
    title(['Measured Velocities per piv ' num2str(i)]);
    colorbar;
    drawnow;
    pause(0.2);
end

figure;
mean_horz_vel = nanmean(horz_vel,3);
mean_vert_vel = nanmean(vert_vel,3);
mean_mag_vel = sqrt(mean_horz_vel.^2 + mean_vert_vel.^2);
cmap = colormap;
cmap(1,:) = [0.85 0.85 0.85];
colormap(cmap)
imagesc(mean_mag_vel, [0 max_vel]); 
colorbar;
title('Mean Speckle Velocities');
xlabel('kernalized pixels in horizontal direction'); ylabel('kernalised pixels in vertical direction');

horz_disp = piv.U;
vert_disp = piv.V;
mean_horz_disp = nanmean(horz_disp,3);
mean_vert_disp = nanmean(vert_disp,3);

%% colormap definition
I = mean_mag_vel;
Ic = round( I/max(max(I))*64);
Ic( Ic == 0) = 1;
Ic(isnan(Ic)) = 1;
colormap hsv;
C = colormap;

%% plotting
log_PDI = mean(abs(BF_filtered),3); 
log_PDI = imresize(log_PDI,[ kgrid.Ny , kgrid.Nx]);% should be done in beamformed dimensions

figure; 

imagesc(20*log10(abs(log_PDI/max(log_PDI(:)))));
colormap gray;
hold on;

lvquiv_mf = 1;
lvquiv_sf = 2;

for indsz = 1:lvquiv_mf:size(piv.X,1);
for indsx = 1:lvquiv_mf:size(piv.X,2);
lvquiv_quiv = quiver(piv.X(indsz,indsx),piv.Z(indsz,indsx),...
    mean_horz_vel(indsz,indsx,1).*lvquiv_sf, mean_vert_vel(indsz,indsx,1).*lvquiv_sf,0,...
    'Color',C(Ic(indsz,indsx),:),'maxheadsize',2);
end
end
hold off;
xlabel('pixels in horizontal direction')
ylabel('pixels in vertical direction ')
title('Alternative Speckle representation')
% figure;
% magn = zeros(kgrid.Ny,kgrid.Nx);
% 
% for i = 1:size(piv.X,1)
%     for j = 1:size(piv.X,2)
%         magn(piv.Z(i,j),piv.X(i,j)) = mean_mag_vel(i,j);
%         
%     end
% end
% 
% imagesc(magn, [0 max_vel]); 
%  colormap jet;
%  colorbar;
% title('Mean Measured Velocities');
% hold on
% for indsz = 1:lvquiv_mf:size(piv.X,1);
% for indsx = 1:lvquiv_mf:size(piv.X,2);
%  quiver(piv.X(indsz,indsx),piv.Z(indsz,indsx),...
%     mean_horz_vel(indsz,indsx,1).*lvquiv_sf, mean_vert_vel(indsz,indsx,1).*lvquiv_sf,0,...
%     'Color',[0 0 0],'maxheadsize',2);
% end
% end
% hold off;

%%

ss1 = piv.ss1s{1};
ss1_vel = ss1(:,:,1,:);
seg_vel = actual_velocities(ss1_vel);
kernels = size(ss1_vel,4);
seg_vel = reshape(seg_vel,[],kernels);
mean_vel = nanmean(seg_vel); % should avoid zero or nan

figure;
cmap = colormap;
cmap(1,:) = [0.85 0.85 0.85];
colormap(cmap)
imagesc(reshape(mean_vel,[size(piv.X,1),size(piv.X,1)]), [0 max_vel]); 
colorbar;

title('Mean Actual Velocities');
xlabel('kernalized pixels in horizontal direction'); ylabel('kernalised pixels in vertical direction');



meas = permute(measured_velocities,[3,1,2]);
meas_vel = reshape(meas,[],kernels);

% Bias = []; Variance = [];
clear Variance Bias
clear mean_meas_vel real_vel
j = 1;
for i = 1:kernels
if ~isnan(mean_mag_vel(i)) && ~isnan(mean_vel(i)) %&& ~isnan(meas_vel(1,i))
Bias(j) =  mean_mag_vel(i) - mean_vel(i) ;
mean_meas_vel(j) = mean_mag_vel(i);
real_vel(j) = mean_vel(i) ;
Variance(j) = nanvar(meas_vel(:,i));
j = j+1;

% Variance = Variance_matrix(:);
end
end
st_dev =[-2.*sqrt(Variance) ; 2.*sqrt(Variance)];



hor_centre = 105;
vert_centre = 120;

clear mean_meas_vel_angle real_vel_angle Variance_angle
angle_inc = 4; i =1;
for theta = 0:angle_inc: 360-angle_inc
    theta_range = [theta , theta + angle_inc -1]; 
    mask_2 = angled_mask(hor_centre,vert_centre,theta_range,kgrid,disc_radius,piv); % wrong use ss1
    bias_angle(i) = nanmean((mean_mag_vel(:) - mean_vel(:)).*mask_2(:));
    mean_meas_vel_angle(i) = nanmean(mean_mag_vel(:).*mask_2(:));
    real_vel_angle(i) = nanmean(mean_vel(:).*mask_2(:));  
    Variance_angle(i) = nanmean(nansum(nanvar(measured_velocities.*mask_2,3)));
    i = i + 1; 
end
Lse_speckle = nanmean((mean_mag_vel(:) - mean_vel(:)).^2);

st_dev_angle =[-2.*sqrt(Variance_angle) ; 2.*sqrt(Variance_angle)];
% debug stats
% h = reshape(horz_vel,[],kernels); v = reshape(vert_vel, [] , kernels);
% mean_velocity = sqrt(nanmean(h(:,i).^2) + nanmean(v(:,i).^2))


% figure;% waarning nan is still there, why so many checks for nans
% plot(mean_meas_vel,'-k');
% hold on 
% plot(st_dev(1,:) + mean_meas_vel,'--k');
% plot(st_dev(2,:) + mean_meas_vel,'--k');
% plot(real_vel,'-r');
% hold off

mean_meas_vel_angle(isnan(mean_meas_vel_angle)) = 0;
st_dev_angle(isnan(st_dev_angle)) = 0;
real_vel_angle(isnan(mean_meas_vel_angle)) = 0; % can be both but needs to be tested

% for clockwise rotation
mean_meas_vel_angle = fliplr(mean_meas_vel_angle);
st_dev_angle = fliplr(st_dev_angle);
real_vel_angle = fliplr(real_vel_angle);

x_axis = angle_inc - 1: angle_inc: 360 - 1;
%% for positive rotation and depth reference

% some elegant way
j = 1;
for  i = 1:length(x_axis)
if x_axis(i) >=90 && x_axis(i) <=270
    clocwise_angles(j) = x_axis(i);
    pstdev(:,j) = st_dev_angle(:,i);
    pmeasvel(j) = mean_meas_vel_angle(i);
    pmeanrealv(j) = real_vel_angle(i);
    j = j + 1;
end
end

figure;% waarning nan is still there, why so many checks for nans

plot(clocwise_angles,pmeasvel,'-k');
hold on 
plot(clocwise_angles,pstdev(1,:) + pmeasvel,'--k');
plot(clocwise_angles,pstdev(2,:) + pmeasvel,'--k');
plot(clocwise_angles,pmeanrealv,'-r');
hold off
xlabel('Clock-wise Angular direction in degrees ');
xlim([90 270]);
ylabel('Velocities(m/s)')
title('Speckle Tracking')
legend('measured','-2*std','+2*std','actual');
% for anti-clockwise rotation
j = 1;
clear anti_clocwise_angles nstdev nmeasvel nmeanrealv
for  i = 1:length(x_axis)
if  x_axis(end - (i-1)) > 270
    anti_clocwise_angles(j) = x_axis(end - (i-1)) - 360;
    nstdev(:,j) = st_dev_angle(:,end - (i-1));
    nmeasvel(j) = mean_meas_vel_angle(end - (i-1));
    nmeanrealv(j) = real_vel_angle(end - (i-1));
    j = j +1;
elseif x_axis(end - (i-1)) <= 90 
    
    anti_clocwise_angles(j) = x_axis(end - (i-1));
    nstdev(:,j) = st_dev_angle(:,end - (i-1));
    nmeasvel(j) = mean_meas_vel_angle(end - (i-1));
    nmeanrealv(j) = real_vel_angle(end - (i-1));
    j = j +1;
end
end

[~,I] = sort(anti_clocwise_angles);

figure;% waarning nan is still there, why so many checks for nans

plot(anti_clocwise_angles(I),nmeasvel(I),'-k');
hold on 
plot(anti_clocwise_angles(I),nstdev(1,I) + nmeasvel(I),'--k');
plot(anti_clocwise_angles(I),nstdev(2,I) + nmeasvel(I),'--k');
plot(anti_clocwise_angles(I),nmeanrealv(I),'-r');
hold off
xlabel('Anti Clock-wise Angular direction in degrees ');
xlim([-90 90]);
ylabel('Velocities(m/s)')
title('Speckle Tracking')
legend('measured','-2*std','+2*std','actual');

%%
% figure;% waarning nan is still there, why so many checks for nans
% 
% plot(x_axis,mean_meas_vel_angle,'-k');
% hold on 
% plot(x_axis,st_dev_angle(1,:) + mean_meas_vel_angle,'--k');
% plot(x_axis,st_dev_angle(2,:) + mean_meas_vel_angle,'--k');
% plot(x_axis,real_vel_angle,'-r');
% hold off
% xlabel('Angular direction in degrees ');
% ylabel('Velcoities(m/s)')
%% per angle visualization later

% direction mean and error

%% Comparison Plots %%
%% Doppler 
dpmask = actual_velocities;
dpmask(dpmask>0)=1;
%dpmask(dpmask==0)=NaN;
figure;imagesc(dpmask);
dopV_grid = imresize(dopV,[kgrid.Ny,kgrid.Nx]);
hor_centre = 105;
vert_centre = 120;

clear mean_meas_vel_angle real_vel_angle Variance_angle
angle_inc = 4; i =1;
for theta = 0:angle_inc: 360-angle_inc
    theta_range = [theta , theta + angle_inc -1]; 
   dmask = Doppler_angled_mask(hor_centre,vert_centre,theta_range,kgrid,disc_radius,piv); % wrong use ss1
%     bias_angle(i) = nanmean((mean_mag_vel(:) - mean_vel(:)).*dmask(:));
    mean_meas_vel_angle(i) = nanmean(dopV_grid(:).*dmask(:).*dpmask(:));
    real_vel_angle(i) = nanmean(actual_velocities(:).*dmask(:));
    i = i+1;
end



% debug stats


mean_meas_vel_angle(isnan(mean_meas_vel_angle)) = 0;

real_vel_angle(isnan(mean_meas_vel_angle)) = 0; % can be both but needs to be tested

% for clockwise rotation
mean_meas_vel_angle = fliplr(mean_meas_vel_angle);

real_vel_angle = fliplr(real_vel_angle);

x_axis = angle_inc - 1: angle_inc: 360 - 1;
%% for positive rotation and depth reference

% some elegant way
j = 1;
for  i = 1:length(x_axis)
if x_axis(i) >=90 && x_axis(i) <=270
    clocwise_angles(j) = x_axis(i);
   
    pmeasvel(j) = mean_meas_vel_angle(i);
    pmeanrealv(j) = -real_vel_angle(i);
    j = j + 1;
end
end

figure;% waarning nan is still there, why so many checks for nans

plot(clocwise_angles,pmeasvel,'-k');
hold on 

plot(clocwise_angles,pmeanrealv,'-r');
hold off
xlabel('Clock-wise Angular direction in degrees ');
xlim([90 270]);
ylabel('Velocities(m/s)')
title('Doppler')
legend('measured','actual');
% for anti-clockwise rotation
j = 1;
clear anti_clocwise_angles nstdev nmeasvel nmeanrealv
for  i = 1:length(x_axis)
if  x_axis(end - (i-1)) > 270
    anti_clocwise_angles(j) = x_axis(end - (i-1)) - 360;
    
    nmeasvel(j) = mean_meas_vel_angle(end - (i-1));
    nmeanrealv(j) = real_vel_angle(end - (i-1));
    j = j +1;
elseif x_axis(end - (i-1)) <= 90 
    
    anti_clocwise_angles(j) = x_axis(end - (i-1));
  
    nmeasvel(j) = mean_meas_vel_angle(end - (i-1));
    nmeanrealv(j) = real_vel_angle(end - (i-1));
    j = j +1;
end
end

[~,I] = sort(anti_clocwise_angles);

figure;% waarning nan is still there, why so many checks for nans

plot(anti_clocwise_angles(I),nmeasvel(I),'-k');
hold on 

plot(anti_clocwise_angles(I),nmeanrealv(I),'-r');
hold off
xlabel('Anti Clock-wise Angular direction in degrees ');
xlim([-90 90]);
ylabel('Velocities(m/s)')
title('Doppler')
legend('measured','actual');