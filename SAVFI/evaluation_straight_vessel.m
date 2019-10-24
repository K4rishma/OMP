function evaluation_straight_vessel(Velocities, n_emissions, savepath)
%-- Example script for evaluating the estimates in a straight vessel tube -
%   converted to function by Jason Voorneveld (j.voorneveld->@->erasmusmc.nl
%

%-- After choosing the specific configuration through acquisition_type,
%-- phantom_type,and flag_display parameters, this script allows computing the
%-- chosen evaluation metrics.

%-- The computed metrics are displayed and performance figures are saved in
%   the folder named "evaluation_results" together with the required *phantom*.mat

%-- Author: Carlos A. Villag√≥mez Hoyos (cavh@elektro.dtu.dk)
%-- $Version: 1.2 $
%-- Added misalignment check - cavh
%-- $Date: 2017/12/12 $
addpath('_aux');

%-- Parameters
if strfind(savepath,'simulation') ~= -1
    acquisition_type = 1;       %-- 1: simulation 
else
    acquisition_type = 2;        %-- 2: measurement
end
if strfind(savepath,'90deg') ~= -1
    phantom_type = 1;           %-- 1: 90deg 
else
    phantom_type = 2;           %-- 2: 105deg
end
flag_save = 1;              %-- 0: do not save plot | 1: save plots


%-- Parse parameter choices
switch acquisition_type
    case 1
        acquisition = 'simulation';
        flag_simu = 1;
    case 2
        acquisition = 'measurement';
        flag_simu = 0;
    otherwise       %-- Otherwise: simulation
        acquisition = 'simulation';
        flag_simu = 1;
end
switch phantom_type
    case 1
        phantom = 'straight_vessel_90deg';
    case 2
        phantom = 'straight_vessel_105deg';
    otherwise       %-- Otherwise: 90deg
        phantom = 'straight_vessel_90deg';
end


%-- Path to the reference grid
reference_path = ['SAVFI',filesep,'reference_grids',filesep,phantom ,filesep, acquisition, filesep,'ref_grid.mat'];


%-- Path to output
[res_dir, res_file, ~] = fileparts(savepath);
path_output_results = fullfile(res_dir, 'evaluation_results');
filename_results=[res_file '_evaluation_results'];


%% Perform evaluation
disp(['Starting evaluation from ',acquisition,' for ',phantom ' and storing them in ' path_output_results])
% Create output directory
if ~exist(path_output_results,'dir')
    if ~mkdir(path_output_results)
        error('Failed to create the output directory!')
    end
end


switch phantom_type
    case 1 	%-- values for evaluating 90 deg vessel
        par.v_peak=0.25;     % Peak velocity [m/s]
        par.z_peak=0.0248;    % Peak position [m] 
        par.d_vessel=0.012;  % Vessel diameter [m]
        par.angle = 90;     % Angle [deg]
        if flag_simu
            par.z_peak=0.025;    % Peak position [m]
            par.angle = 90;   % Angle [deg] (Opposite direction)
        end
        
    case 2 	%--  values for evaluating 105 deg vessel
        par.v_peak=0.25;               % Peak velocity [m/s]
        par.z_peak=0.0248;              % Peak position [m] 
        par.d_vessel=0.012/cosd(15);   % Vessel diameter [m]
        par.angle = 75;               % Angle [deg]
        if flag_simu
            par.z_peak=0.025;              % Peak position [m]
            par.angle = 105;            % Angle [deg]  % Opposite direction
        end
        
    otherwise       %-- Do deal with bad values
        par.v_peak=0.25;     % Peak velocity [m/s]
        par.z_peak=0.025;    % Peak position [m] 
        par.d_vessel=0.012;  % Vessel diameter [m]
        par.angle = 90;      % Angle [deg]
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
v_true=(reference.flow_grid.z-(par.z_peak-0.5*par.d_vessel)).*(reference.flow_grid.z-(par.z_peak+0.5*par.d_vessel));
v_true(reference.flow_grid.z<(par.z_peak-0.5*par.d_vessel))=0;
v_true(reference.flow_grid.z>(par.z_peak+0.5*par.d_vessel))=0;
v_true=par.v_peak*v_true./min(v_true);
d_true=ones(size(v_true))*par.angle;

%-- Convert data to polar (magnitude and angle)
Velocities.magnitude=sqrt(Velocities.x.^2+Velocities.y.^2+Velocities.z.^2);
Velocities.phi=atan2d(Velocities.y,Velocities.x);
Velocities.theta=atan2d(Velocities.x,Velocities.z);  % changed here - was not correct for my use case


%-- Verify alignment by cross-correlation
coeff=squeeze(xcorr(mean(Velocities.magnitude,4),v_true)); % Xcorr the true profile with mean profile
coeff_center_idx=(size(coeff, 1)-1)/2+1;
[bla idx_center] = max(coeff);
center_displacement_samples=idx_center-coeff_center_idx;
if (center_displacement_samples~=0)
    % If more than one sample misaligned warn the participant that samples
    % are not aligned, and estimate a compensation shift that he can use
    % when estimating the velocity grid.
    
    % Coarse distance
    coarse_adjustment=center_displacement_samples*(reference.flow_grid.z(2)-reference.flow_grid.z(1));
    % Parabolic interpolation for fine adjustment
    fine_adjustment= (reference.flow_grid.z(2)-reference.flow_grid.z(1))*(coeff(idx_center+1 )-coeff(idx_center-1))...
        ./(2*(coeff(idx_center+1)-2*coeff(idx_center)+coeff(idx_center-1)));
    % Compensation adjustement
    final_adjustment= coarse_adjustment+fine_adjustment;
    warning(['Important: The results seem misaligned, this can affect the evaluation. Please adjust estimation grid by ' num2str(final_adjustment*1000) 'mm'] )
    
end

%-- Generate statistics
vmag_est_mean=squeeze(mean(Velocities.magnitude, 4))';
vmag_est_std=squeeze(std(Velocities.magnitude, [], 4))';

%-- Circular statistics for angle
r=sum(exp(1i*deg2rad(Velocities.theta)),4);   % Map to unit circle and sum
angle_est_mean= squeeze(rad2deg(angle(r)))';  % The angle is the mean
angle_est_mean(angle_est_mean<0)= ...
    angle_est_mean(angle_est_mean<0)+360;     % make it 0 - 360 range

%-- angular deviation according to Zar  (26.20). (Biostatistical Analysis, J. H. Zar)
angle_est_std = squeeze(rad2deg(real(sqrt(2*(1-(abs(r)/size(Velocities.theta,4)))))))';

%-- Estimate the error
v_err=vmag_est_mean-v_true;
angle_err=angle_est_mean-d_true;

%-- Make NaN everything outside the vessel boundaries
pct2consider=0.5;
v_err(reference.flow_grid.z<(par.z_peak-pct2consider*par.d_vessel))=NaN;
v_err(reference.flow_grid.z>(par.z_peak+pct2consider*par.d_vessel))=NaN;
angle_err(reference.flow_grid.z<(par.z_peak-pct2consider*par.d_vessel))=NaN;
angle_err(reference.flow_grid.z>(par.z_peak+pct2consider*par.d_vessel))=NaN;
vmag_est_mean(reference.flow_grid.z<(par.z_peak-pct2consider*par.d_vessel))=NaN;
vmag_est_mean(reference.flow_grid.z>(par.z_peak+pct2consider*par.d_vessel))=NaN;
vmag_est_std(reference.flow_grid.z<(par.z_peak-pct2consider*par.d_vessel))=NaN;
vmag_est_std(reference.flow_grid.z>(par.z_peak+pct2consider*par.d_vessel))=NaN;
angle_est_mean(reference.flow_grid.z<(par.z_peak-pct2consider*par.d_vessel))=NaN;
angle_est_mean(reference.flow_grid.z>(par.z_peak+pct2consider*par.d_vessel))=NaN;
angle_est_std(reference.flow_grid.z<(par.z_peak-pct2consider*par.d_vessel))=NaN;
angle_est_std(reference.flow_grid.z>(par.z_peak+pct2consider*par.d_vessel))=NaN;

%-- Make the statitistics relative to the peak value in percentage
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
set(0,'defaultaxesfontsize',15)
set(0,'defaulttextfontsize',15);
set(0,'defaultlinelinewidth',2)
figure('position', [0, 0, 800, 800])
subplot(3,1,[1 2])
plot(reference.flow_grid.z*1000, v_true, 'r')
hold on;
shadedErrorBar(reference.flow_grid.z*1000,vmag_est_mean,vmag_est_std,'k',1);
hold off;
ylabel('Estimated velocity [m/s]')
title('Velocity magnitude')
grid on
subplot(3,1,3)
plot(reference.flow_grid.z*1000, d_true, 'r')
hold on;
shadedErrorBar(reference.flow_grid.z*1000,angle_est_mean,angle_est_std,'k',1);
xlabel('Axial position [mm]')
ylabel('Estimated angle [deg]')
grid on
% display scores on image

velmagtext = {'Vel. mag. %:',velmagstr};
annotation('textbox',[0.7 0.7 0.18 0.2],'EdgeColor','None', 'String',velmagtext);
velangtext = ['Vel. ang. ∞:', velangstr];
annotation('textbox',[0.15 0.12 0.75 0.2],'EdgeColor','None', 'String',velangtext) 

% Save figure if enabled
if flag_save
    path_output_log=[path_output_results filesep filename_results];
    disp(['Image saved in "',path_output_log, '.png','"'])
    print(path_output_log ,'-dpng')
    disp(['Results saved in "',path_output_log, '.mat','"'])
    results = struct();
    results.Velocities = Velocities;
    results.n_emissions = n_emissions;
    save(path_output_log,'-struct','results');
end

disp('Evaluation Done')
