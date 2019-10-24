close all hidden
clear all
clc

%% Load the data
dir_name = 'K:\Data\Pieter\';

% Load the parameter file
% load([dir_name '\Parameters.mat'])
% load([dir_name '\BF_5000.mat'])

Ne_new = 25;
load('C:\Users\k.k\Desktop\wetransfer-e90f62\spectrum_analysis\BeamformedDataForKarishma.mat');
BF = BF(:,:,1:Ne_new); 
% load dat file
%fid = fopen([dir_name '\BF'], 'rb');
%bfr = fread(fid, [1 2*BFrecon.Nx*BFrecon.Nz*Ne_new], 'single=>single');
%fclose(fid);
%BF = complex(bfr(1:2:end),bfr(2:2:end));
%BF = reshape(BF,BFrecon.Nz,BFrecon.Nx, Ne_new);

Nz_start = 30;                  % Crop the image starting from Nz_start
Nz_end = 209;                  % Crop the image starting from Nz_start
Nx_start = 1;
Nx_end = 128;
BF = BF(Nz_start:Nz_end,Nx_start:Nx_end,:);    % Crop the dataset
Nz = size(BF,1);         % Update Nz
Nx = size(BF,2);         % Update Nx


%% Doppler Filtering
% FsDoppler = param.dopPRF;
% f1 = 40;
% filter_order = 8;
% [b,a] = butter(filter_order,f1/(FsDoppler/2),'high');
% 
% BF_filtered = single(filtfilt(b,a,double(reshape(BF,Nz*Nx,[])')))';
% BF_filtered = reshape(BF_filtered,Nz,Nx,Ne_new);

% [BF_filtered,U, Sigma, V, lambdas1, lambdas2] = svdfilterAuto(BF);
[BF_filtered,~, ~, ~] = svdfilter(BF, 7, 0);


BFN = reshape(BF,Nz*Nx,[]);
% power = abs(BFN(1:50,1:150)*BFN(1:50,1:150)');
% figure;imagesc(power);
% BF_filtered = BF_filtered.^8;



% BFN = reshape(BF,Nz*Nx,[]);
% [U,S,V] = svdecon(BFN);
% S(:,1:250) = 0;
% BFN = U*S*V';
% BFN = reshape(BFN,size(BF));
%% Do the oversampling
lambda = 60e-6;
dx_old = median(diff(BFrecon.x_axis)).*1e-3;
dz_old = median(diff(BFrecon.z_axis)).*1e-3;

dx_new = lambda / 2 / 2;
dz_new = dx_new;         % Square pixel size 

Nx_new = round( Nx * (dx_old / dx_new));
Nx_new = Nx_new + rem(Nx_new,2);
Nz_new = round( Nz * (dz_old / dz_new));
Nz_new = Nz_new + rem(Nz_new,2);

BF_new = ifft2(fftshift(fft2(BF_filtered,Nz,Nx),2),Nz_new,Nx_new);

figure; imagesc(abs(mean(BF_new,3)));
%% autocorrelation

auto1 = sum(conj(BF_new(:,:,1:end-1)).*(BF_new(:,:,2:end)),3)/(size(BF_new,3)-1);
dopV = (param.c * param.dopPRF * angle(auto1))/(4*BFrecon.FsIQ*pi);

figure;imagesc(dopV)


%% Do the PIV
scan = ImageProperties();

scan.frames = permute(BF_new, [1,2,4,3]);
scan.dx = dx_new*1e3;  %mm
scan.dz = dz_new*1e3;  %mm
scan.origin = [BFrecon.z_axis(Nz_start), BFrecon.x_axis(Nx_start)];
scan.FR = param.dopPRF;
scan.update()

scan.create_bmode_by_coherent_compounding();
scan.ensemble_average_bmode();
%scan.play_frames

temp = squeeze(mean(abs(scan.frames),4));  % get PWD
% temp = temp - max(temp(:)); %norm to max

scan.imMask_cart = temp > 40;
%figure; imshow(scan.imMask_cart);
scan.ROI = [1, 1, scan.frame_size(2), scan.frame_size(1)];
scan.play_bmode_with_mask


% scan.play_frames

piv = PIV();
    
fdims = [size(scan.frames(:,:,1,1)), Opts.n_ave*scan.num_ang+Opts.nm_ave-1];
piv.init(scan.num_frames, fdims);

figure;
for i=1:1:Ne_new
    I = imagesc(abs(scan.frames(:,:,1,i)),[0 300]);%,[0 0.005]); 
%     I = imadjust(abs(scan.frames(:,:,1,i)), [0.3 1], [0 1]);
%     imshow(I);
    title(['Beamformed image ', num2str(i)]);
    colormap gray
    colorbar;drawnow;
     pause(0.4);
end

video = VideoWriter('newfile.avi', 'Uncompressed AVI');
video.FrameRate = 2;
open(video)
for k=1:Ne_new
    
   W = abs(scan.frames(:,:,1,k));
   %W=reshape(w,y,x); %Where W is a Matrix
   lsd=imagesc(W);
   title(['Frame ', num2str(k)]);
   colormap gray;
   if k == 1;
       C = colorbar;
   else
       colorbar = C;
   end
   drawnow
   F(k) = getframe(gcf);
   writeVideo(video,F(k))
end
close(video)


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

pivfilt = PostProcess(piv, scan);
pivfilt.postprocess();
pivfilt.convert_disp_to_vel();
pivfilt.play_velmag();
scan.PIV_view_cart(pivfilt, [-60 0],1.0);

%% Visualisation
close all
PDI = mean(abs(BF_new),3);

z_axis_show = linspace(BFrecon.z_axis(Nz_start),BFrecon.z_axis(Nz_end),Nz_new);
x_axis_show = linspace(BFrecon.x_axis(Nx_start),BFrecon.x_axis(Nx_end),Nx_new);

figure
subplot(121)
pdimax = max(PDI(:));
PDIshow = 20*log10(PDI./pdimax);
imagesc(x_axis_show,z_axis_show,PDIshow,[-30 0])
xlabel('Width [mm]')
ylabel('Depth [mm]')
title('Power Doppler Image')
colormap(hot);
colorbar
axis equal
axis tight

subplot(122)
for i = 1:Ne_new
    IMshow = abs(BF_new(:,:,i));
    imagesc(x_axis_show,z_axis_show,IMshow,[0 100])
    xlabel('Width [mm]')
    ylabel('Depth [mm]')
    title('Power Doppler Image')
    colormap(hot);
    colorbar
    axis equal
    axis tight
    drawnow;
end

