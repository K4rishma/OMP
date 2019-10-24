function evaluation_spinning_disk(Velocities, n_emissions, savepath)

%-- Example script to evaluate the estimates in a spinning disk
%   converted to function by Jason Voorneveld (j.voorneveld->@->erasmusmc.nl

%-- By default, the computed metrics are saved in a text file in the folder named "evaluation"

%-- Author: Carlos A. Villagómez Hoyos (cavh@elektro.dtu.dk)
%-- $Version: 1.0 $
%-- $Date: 2017/11/14 $
addpath('_aux');

%-- Parameters
flag_save = 1;           %-- 0: do not save plot | 1: save plots

phantom = 'spinning_disk';
acquisition = 'simulation';

%-- Phantom parameters
par.center=[0 0 25]/1000;       % center [m]
par.v_peak=[0.25];              % peak velocity [m/s]
par.radius=0.0075;              % radius [m]

%-- Evaluation parameters
par.mask_ini=0.1;                 % Mask initial radius in percent
par.mask_end=1;                   % Mask final radius in percent

%-- Path to the reference grid
reference_path = ['SAVFI',filesep,'reference_grids',filesep,phantom ,filesep, acquisition, filesep,'ref_grid.mat'];

%-- Path to output
[res_dir, res_file, ~] = fileparts(savepath);
path_output_results = fullfile(res_dir, 'evaluation_results');
filename_results=[res_file '_evaluation_results'];


%%  Perform evaluation
disp(['Starting evaluation for ',phantom])
% Create output directory
if ~exist(path_output_results,'dir')
    if ~mkdir(path_output_results)
        error('Failed to create the output directory!')
    end
end

%- Load the reference grid
reference=load(reference_path);

%- make it an option to not process the entire dataset
if size(Velocities.x,4) < length(reference.flow_grid.t)
    reference.flow_grid.t = reference.flow_grid.t(1:size(Velocities.x,4));
end

%- Verify that the results matches the reference flow_grid
grid_names=['x' 'y' 'z' 't'];
ref_size = [length(reference.flow_grid.(grid_names(1))) ...
    length(reference.flow_grid.(grid_names(2))) ...
    length(reference.flow_grid.(grid_names(3))) ...
    length(reference.flow_grid.(grid_names(4)))];
for it_grid=1:length(grid_names)-1
    if ~issame(ref_size,size(Velocities.(grid_names(it_grid))))
        error(['The data dimension "' grid_names(it_grid) '" does not match the reference flow_grid.' ]);
    end
end

if exist('n_emissions','var') && isa(n_emissions,'numeric')
    par.n_emissions=real(abs(n_emissions));
else
    error('The number of emissions is needed for estimating the results.')
end


%-- Generate true velocities magnitude and direction
[X Z] = meshgrid(reference.flow_grid.x,reference.flow_grid.z);
velocity_true= par.v_peak*(sqrt((X-par.center(1)).^2+(Z-par.center(3)).^2)./par.radius);
angles_true=-1*(radtodeg(atan2(-1*(Z-par.center(3)),(X-par.center(1)))))+180;
mask_idx=find((velocity_true>=(par.v_peak*par.mask_end)) + (velocity_true<=(par.v_peak*par.mask_ini))  );
angles_true(mask_idx)=NaN;
velocity_true(mask_idx)=NaN;

%-- Convert data to polar (magnitude and angle)
Velocities.magnitude=sqrt(Velocities.x.^2+Velocities.y.^2+Velocities.z.^2);
Velocities.theta=atan2d(-1*Velocities.x,Velocities.z)+180;

%-- Generate statistics
vmag_est_mean=squeeze(mean(Velocities.magnitude, 4))';
vmag_est_std=squeeze(std(Velocities.magnitude, [], 4))';

%-- Circular statistics
r=sum(exp(1i*deg2rad(Velocities.theta)),4);   % Map to unit circle and sum
angle_est_mean= squeeze(rad2deg(angle(r)))';  % The angle is the mean
angle_est_mean(angle_est_mean<0)= ...
    angle_est_mean(angle_est_mean<0)+360;     % make it 0 - 360 range

%-- angular deviation according to Zar  (26.20). (Biostatistical Analysis, J. H. Zar)
angle_est_std = real(squeeze(rad2deg(sqrt(2*(1-(abs(r)/size(Velocities.theta,4)))))))';


%-- Estimate the error
v_err=abs(vmag_est_mean-velocity_true);
angle_err=exp(1i*deg2rad(angle_est_mean))+exp(-1i*deg2rad(angles_true));  % Map to unit circle and substract phases
angle_err=abs(real(angle_err))+1i*(imag(angle_err));  % Make it two quadrants only
angle_err= abs(rad2deg(2*angle(angle_err)));   % Multiply by 2 to make it 360


%-- Apply mask to estimates
vmag_est_mean(mask_idx)=NaN;
vmag_est_std(mask_idx)=NaN;
angle_est_mean(mask_idx)=NaN;
angle_est_std(mask_idx)=NaN;

%-- Make the statitistics relative to the peak value for percentage
v_bias=nanmean(v_err(:))/par.v_peak*100;
v_std=sqrt(nanmean(vmag_est_std(:).^2))/par.v_peak*100;
weight_v_std=v_std*sqrt(1+(par.n_emissions/5));
ang_bias=nanmean(angle_err(:));
ang_std=sqrt(nanmean(angle_est_std(:).^2));
weight_ang_std= ang_std*sqrt(1+(par.n_emissions/5));


%-- Display scores
fprintf('--------------------------------------------------\n');
fprintf('Velocity magnitude metrics [pct respect %f m/s, no_emissions = %f ]: \n',par.v_peak,par.n_emissions);
velmagstr = sprintf('Bias: %5.3f, std: %5.3f wstd: %5.3f \n', v_bias, v_std, weight_v_std);
fprintf(velmagstr);
fprintf('Velocity angle metrics [deg, no_emissions = %f ]: \n',par.n_emissions);
velangstr = sprintf('Bias: %5.3f, std: %5.3f wstd: %5.3f \n', ang_bias, ang_std, weight_ang_std);
fprintf(velangstr);
fprintf('--------------------------------------------------\n');






%% Display figures
set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultfigureposition',[0, 0, 700, 600]);   % Set aspect ratio
path_output_log=[path_output_results filesep filename_results];

% First figure. The magnitude of the velocity
figure();
contourf(reference.flow_grid.x*1000,reference.flow_grid.z*1000, vmag_est_mean, [0:0.025:1.1]*par.v_peak,':','LineWidth',1);
colormap(jet(length(0:0.025:1.1)));
cb=colorbar;
hold on;
ylabel(cb,'Magnitude [m/s]')
[C,h] = contour(reference.flow_grid.x*1000,reference.flow_grid.z*1000,velocity_true,[0:0.2:1]*par.v_peak,'w--','LineWidth', 3);
clabel(C,h,[0.1:0.1:1]*par.v_peak,'FontSize',18,'Color','white')
title(['Estimated velocity [m/s]'])
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
axis xy
axis equal
set(gca,'YDir','Reverse')
if flag_save
    disp(['Velocity magnitude image saved in "',[path_output_log '_Vmag'], '.png','"'])
    print([path_output_log '_Vmag'] ,'-dpng')
end


% Second figure. The velocity error
figure;
contourf(reference.flow_grid.x*1000, reference.flow_grid.z*1000, v_err/par.v_peak*100, [0:2:20]);
colormap(thermal(10))
c = colorbar;
ylabel(c,'Relative error [%]')
title(['Velocity magnitude error [%]'])
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
axis xy
axis equal
set(gca,'YDir','Reverse')
if flag_save
    disp(['Velocity error image saved in "',[path_output_log '_Verror'], '.png','"'])
    print([path_output_log '_Verror'] ,'-dpng')
end

% Third figure. The velocity standard deviation
figure;
contourf(reference.flow_grid.x*1000, reference.flow_grid.z*1000, vmag_est_std/par.v_peak*100, [0:2:20]);
colormap(thermal(10))
c = colorbar;
ylabel(c,'Relative Std. Dev [%]')
title(['Velocity Std. Dev [%]'])
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
axis xy
axis equal
set(gca,'YDir','Reverse')
if flag_save
    disp(['Velocity standard deviation image saved in "',[path_output_log '_Vstd'], '.png','"'])
    print([path_output_log '_Vstd'] ,'-dpng')
end

% Fourth figure.  The angles mean
figure();
pcolor(reference.flow_grid.x*1000, reference.flow_grid.z*1000, angle_est_mean);
colormap(hsv(360))
shading flat
c = colorbar;
ylabel(c,'Angle mean [\circ]')
hold on;
contour(reference.flow_grid.x*1000, reference.flow_grid.z*1000, angle_est_mean, [0:15:360],'k:','LineWidth',1);
[C,h] = contour(reference.flow_grid.x*1000, reference.flow_grid.z*1000,angles_true,[30 60 90 120 150 180 210 240 270 300 330],'w--','LineWidth', 2);
clabel(C,h,[ 30 60 90 120 150 180 210 240 270 300 330],'Color','white','FontSize',18,'LabelSpacing',296)
title(['Estimated Angle [\circ]'])
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
axis xy
axis equal
set(gca,'YDir','Reverse')
if flag_save
    disp(['Velocity angle image saved in "',[path_output_log '_Angle'], '.png','"'])
    print([path_output_log '_Angle'] ,'-dpng')
end

% Fifth figure. The angle error
figure;
contourf(reference.flow_grid.x*1000, reference.flow_grid.z*1000, angle_err, [0:2:20]);
colormap(thermal(10))
c = colorbar;
ylabel(c,'Angle error [\circ]')
title(['Angle error [\circ]'])
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
axis xy
axis equal
set(gca,'YDir','Reverse')
if flag_save
    disp(['Velocity angle error image saved in "',[path_output_log '_Aerror'], '.png','"'])
    print([path_output_log '_Aerror'] ,'-dpng')
end


% Sixth figure. The angle standard deviation
figure;
contourf(reference.flow_grid.x*1000, reference.flow_grid.z*1000, angle_est_std , [0:2:20]);
colormap(thermal(10))
c = colorbar;
ylabel(c,'Angle Std. Dev [\circ]')
title(['Angle Std. Dev[\circ]'])
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
axis xy
axis equal
set(gca,'YDir','Reverse')

if flag_save
    disp(['Velocity angle standard deviation image saved in "',[path_output_log '_Astd'], '.png','"'])
    print([path_output_log '_Astd'] ,'-dpng')
end

if flag_save
   
    disp(['Results saved in "',path_output_log, '.mat','"'])
    results = struct();
    results.Velocities = Velocities;
    results.n_emissions = n_emissions;
    save(path_output_log,'-struct','results');
end


disp('Evaluation Done')




