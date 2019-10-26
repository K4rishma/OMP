% filtering removed % temp mask(x,rot) % change
% angleFromDerivatives %loadfile
close all hidden
clear all
clc
 addpath('statsLib');  % this assumes 2DEchoPIV folder is working dir.
 addpath('Plotting');  % this assumes 2DEchoPIV folder is working dir.
 addpath('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\k-wave-toolbox-version-1.2.1\k-Wave');
%% Load the data
 close all;
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTracking -10 0 10 Ne 150.mat')
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTracking0_Ne_20_4e5_dx_.mat')
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTracking-10  -5   0   5  10_Ne_50_4e5_dx_.mat')
%  load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\parameters_angles_speed2.2096ms_sametimecompound_slowScan_degree_-10   0  10_Ne_100.mat')


% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\Circularsmallgridzerokerfhighsensors-10   0  10_Ne_30_4e5_dx_.mat')

%% working
%  load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTrackingHighVelocity-10  -5   0   5  10_Ne_35_4e5_dx_.mat')

%       load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTrackingHighVelocityVessel-10  -5   0   5  10_Ne_35_4e5_dx_.mat')
%    load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTrackingHighVelocityDisc-10  -5   0   5  10_Ne_35_4e5_dx_.mat')
%   load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTrackingHighVelocityHalfDisc-10  -5   0   5  10_Ne_35_4e5_dx_.mat')

%     load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\SpeckleTrackingHighVelocityDiscdensitylowerConv0_Ne_25_4e5_dx_.mat')
%    load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\twicethinDiscRBC0_Ne_100_4e5_dx_.mat')
%  load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\twicethinDiscRBC0_Ne_30_4e5_dx_.mat');
%    velo =  1m/2 when vd = 2.8m/s(true vel is not 1m/s))
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\twicethinDiscRBCprf0_Ne_30_4e5_dx_.mat')
% Nt increase to decrease the max doppler velocity
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\twicethinDiscRBCprfconv0_Ne_30_4e5_dx_.mat')
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\twicethinDiscRBCprfconvwt0_Ne_60_4e5_dx_.mat')
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\twicethinDiscRBCprfconvwtcentre0_Ne_30_4e5_dx_.mat')
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\twicethinDiscRBCprfconvwtcentre-10  -5   0   5  10_Ne_30_4e5_dx_.mat')
% good(limit is set) result for mean vertical and horizontal velocities, high mean
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\straight_vessel-10  -5   0   5  10_Ne_30_4e5_dx_.mat')

% reduce the velocity and start straight vessels,,,,,
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\straight_vessel_freqRemov_lowresol0_Ne_25_4e5_dx_.mat')
% low resoution on no smoothening of medium speed and density
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\straight_vessel_freqRemov_lowresol-10  -5   0   5  10_Ne_25_4e5_dx_.mat')
% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\straight_vessel_freqRemov_lowresol_densityg-10  -5   0   5  10_Ne_30_4e5_dx_.mat')

% load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\straight_vessel_freqRemov_lowresol_densityg-10  -5   0   5  10_Ne_60_4e5_dx_.mat')
 load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\Slantedmin0.3-10  -5   0   5  10_Ne_60_4e5_dx_.mat')
%  load('C:\Users\k.k\Desktop\wetransfer-e90f62\k_wave_analysis\datasets\cylindrical_H\horizontal0.3-10  -5   0   5  10_Ne_60_4e5_dx_.mat')
% observations with 1 as start scan, filtering operation is evident in 10, 
 
IM = (data{1,6});
%IM = Rfstore;
kgrid = data{1,4};
Scan = data{1,1};
Trans = data{1,2};
medium = data{1,3};
CI = data{1,7};

Ne_new = size(data{1,7},3);

% for comparing it with Pieter's code
% BF = CI; % for another script
% BF = BF(:,:,1:Ne_new); 
BF1 = data{1,6};
 clear BF
 BF(:,:,1:Scan.Ne) = sum(BF1(:,:,:,1:Scan.Ne),3);


Nz_start = 1;                  % Crop the image starting from Nz_start
Nz_end = size(BF,1);                  % Crop the image starting from Nz_start
Nx_start = 1;
Nx_end = size(BF,2);
BF = BF(Nz_start:Nz_end,Nx_start:Nx_end,1:Ne_new);    % Crop the dataset
Nz = size(BF,1);         % Update Nz
Nx = size(BF,2);         % Update Nx

%% ground truth - velocities
SimTime = kgrid.Nt*kgrid.dt;
c0 = medium.sound_speed(1,1);
velocity = c0/(SimTime*Trans.f0*4)
%% Doppler Filtering
FsDoppler = 1/(kgrid.Nt*kgrid.dt);
f1 = 0.5/30*FsDoppler;
filter_order = 8;
[b,a] = butter(8,f1/(FsDoppler/2),'high');
BF_filtered = single(filtfilt(b,a,double(reshape(BF,Nz*Nx,[])')))';
BF_filtered = reshape(BF_filtered,Nz,Nx,Ne_new);
 BF_filtered = BF;
% BF_new = BF_filtered;
% [BF_filtered,~, ~, ~] = svdfilter(BF,8, 0);
% 
% BFN = reshape(BF,Nz*Nx,[]);
% power = abs(BFN(1:50,1:150)*BFN(1:50,1:150)');
% figure;imagesc(power);

%% Do the oversampling
% lambda = 192.5e-6;
% dx_old = kgrid.dx;
% dz_old = kgrid.dy;
% 
% dx_new = lambda / 2 / 2;
% dz_new = dx_new;        
% 
% Nx_new = round( Nx * (dx_old / dx_new));
% Nx_new = Nx_new + rem(Nx_new,2);
% Nz_new = round( Nz * (dz_old / dz_new));
% Nz_new = Nz_new + rem(Nz_new,2);
% 
% BF_new = ifft2(fftshift(fft2(BF_filtered,Nz,Nx),2),Nz_new,Nx_new);
% % BF_new = imresize(BF,[Nz_new,Nx_new]);

figure; imagesc(mean(abs(BF_filtered),3));title('filtered PDI');
%% autocorrelation
k = 10;
auto1 = sum(conj(BF_filtered(:,:,1:end-1)).*(BF_filtered(:,:,2:end)),3)./(size(BF_filtered,3)-1);
%dopV = (param.c * param.dopPRF * angle(auto1))/(4*BFrecon.FsIQ*pi);
dopV = (medium.sound_speed(1,1)/(kgrid.Nt*kgrid.dt).* angle(auto1))./(4*Trans.f0*pi);
figure;imagesc(dopV);title('doppler image');

%put the ustb and visualisation of pcc
%% Do the PIV
scan = ImageProperties();
clear BF_new
for j = 1:Ne_new
BF_new(:,:,j) = imresize(BF_filtered(:,:,j),[kgrid.Nx kgrid.Ny]);
end
% BF_new = imresize(BF_new, [kgrid.Ny kgrid.Nx Scan.Ne]);
scan.frames = permute(BF_new, [1,2,4,3]);% IM is not useful as processing is not based on angles
scan.dx = kgrid.dy*1e3;%dx_new*1e3;  %mm
scan.dz = kgrid.dx*1e3;%dz_new*1e3;  %mm
scan.origin = [0, 0 ];
scan.FR = 1/(kgrid.Nt*kgrid.dt);
scan.update()

scan.create_bmode_by_coherent_compounding();
scan.ensemble_average_bmode();
%scan.play_frames

temp = squeeze(mean(abs(scan.frames),4))./nanmax(nanmax(squeeze(mean(abs(scan.frames),4))));%,'all');  % get PWD
% temp = temp - max(temp(:)); %norm to max


% temp = temp > graythresh(temp);% + 0.60*graythresh(temp);
% figure; imagesc(temp); title('temp');
% temp(:,1:81) = 0;
% temp(:,119:end) = 0;

% V1_x = kgrid.Nx/2 - 30: 160;     
 V1_x = kgrid.Nx/2 - 50: 135;%kgrid.Nx/2 + 50; 150     % for slanted vessel
% V1_x = kgrid.Nx/2 - 50: 150;%kgrid.Nx/2 + 50; 
V1_y = kgrid.Ny/2 - 6: 1: kgrid.Ny/2 + 6;
temp = 0.*temp;
temp(V1_x, V1_y) = 1;
temp = imrotate(temp,-60,'crop');
figure; imagesc(temp); title('temp');
scan.imMask_cart = temp;
figure; imshow(scan.imMask_cart);
dopV_grid = imresize(dopV,[kgrid.Ny,kgrid.Nx]);
figure; imagesc(dopV_grid.*scan.imMask_cart);
scan.ROI = [1, 1, scan.frame_size(2), scan.frame_size(1)];
scan.play_bmode_with_mask


% scan.play_frames

piv = PIV();
    
fdims = [size(scan.frames(:,:,1,1)), Opts.n_ave*scan.num_ang+Opts.nm_ave-1];
piv.init(scan.num_frames, fdims);



% Process PIV - TODO - I hate that this loop is not in the class
tic;
h = waitbar(0,'Processing PIV');
for i = piv.piv_start:piv.piv_stop
    piv.image1 = scan.get_ref_image_stack(i);
    piv.image2 = scan.get_target_image_stack(i);

    piv.process (scan, i);

    if Opts.corr_over_quiv
        piv.live_view_corr_update(); 
    end
    if Opts.live_view
        piv.live_view_quiv_update();
    end
    waitbar((i-piv.piv_start)/(piv.piv_stop-piv.piv_start),h);
end            
timer = toc;
close(h);
piv.clean_up();
fprintf('Time to process PIV = %4.2f s or %4.2f fps \n' ,timer, (piv.piv_stop-piv.piv_start)/timer)


ind1 = piv.ss1s;
ind2 = piv.ss2s;

i1 = squeeze(ind1(:,:,1,:));
i2 = squeeze(ind2(:,:,1,:));

% rectangle('Position',[i1(),i1,5,10],'FaceColor',[0 .5 .5],'EdgeColor','b',...
%     'LineWidth',3)

pivfilt = PostProcess(piv, scan);
%pivfilt.postprocess();
pivfilt.convert_disp_to_vel();
pivfilt.play_velmag();
scan.PIV_view_cart(pivfilt, [-60 0],1.0);

% Visualisation
%close all
% PDI = mean(abs(BF_new),3);
% 
% 
% figure
% %subplot(121)
% pdimax = max(PDI(:));
% PDIshow = 20*log10(PDI./pdimax);
% imagesc(PDIshow)
% xlabel('Width [mm]')
% ylabel('Depth [mm]')
% title('Power Doppler Image')
% colormap(hot);
% colorbar
% axis equal
% axis tight
% 
% BF_new = BF_filtered;
% figure;
% for i = 1:1:piv.piv_stop;%Scan.Ne
%     IMshow = mean(real((BF_new(:,:,Opts.n_ave*(i-1)+1:...
%                 i*Opts.n_ave))),3);
% %      IMshow = abs(real(BF_new(:,:,i)));
%     imagesc(IMshow)
%     xlabel('Width [mm]')
%     ylabel('Depth [mm]')
%     title('Power Doppler Image')
%     colormap(hot);
%     colorbar
%     axis equal
%     axis tight
%     drawnow;
%     pause(1)
%     
% end

% imtool(IMshow);
s1 = piv.ss1s; s2 = piv.ss2s;
p = cell2mat(s1);
q1=p(:,:,1,1);
p1 = cell2mat(s2);
q2=p1(:,:,1,1);
