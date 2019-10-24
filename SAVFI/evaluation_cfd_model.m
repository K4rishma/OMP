function evaluation_cfd_model(Velocities, savepath)

%-- Example script to evaluate the estimates from a singe time instance in a CFD model.
%  converted to function by Jason Voorneveld (j.voorneveld->@->erasmusmc.nl

% The CFD model used for this example is publicly available at
% http://www.biommeda.ugent.be/biomedical-ultrasound-research/download-datasets-scatterer-phantoms
% and the model is fully described in:
% - Ultrasound simulation of complex flow velocity fields based on computational fluid dynamics.
%   A. Swillens, L. Løvstakken, J. Kips, H. Torp, and P. Segers.
%   IEEE Trans. Ultrason., Ferroelec., Freq. Contr., 56(3):546–556, 2009.

%-- By default, the computed metrics are saved in a text file in the folder named "evaluation"

%-- Author: Carlos A. Villagómez Hoyos (cavh@elektro.dtu.dk)
%-- $Version: 1.0 $
%-- $Date: 2017/11/14 $
addpath('_aux');

%-- Parameters
flag_save = 1;              %-- 0: do not save plot | 1: save plots

phantom = 'carotid_bifurcation';
acquisition = 'simulation';
if contains(savepath,'testing')
    datatype = 'testing';
else
    datatype = 'training';
end

%-- Path to the reference grid
reference_path = ['SAVFI',filesep,'reference_grids',filesep,phantom ,filesep, datatype, filesep, acquisition, filesep,'ref_grid.mat'];

%-- Path to output
[res_dir, res_file, ~] = fileparts(savepath);
path_output_results = fullfile(res_dir, 'evaluation_results');
filename_results=[res_file '_evaluation_results'];



%% Perform evaluation
disp(['Starting evaluation for ',phantom])
% Create output directory
if ~exist(path_output_results,'dir')
    if ~mkdir(path_output_results)
        error('Failed to create the output directory!')
    end
end

%-- Phantom reference velocities
ref_velocities_path=['_aux' filesep 'ref_CFD_velocities.mat'];

%-- Load reference velocities
reference=load(ref_velocities_path);

%- Load the reference grid
temp=load(reference_path);
reference.flow_grid=temp.flow_grid;


%- Verify that the results matches the reference flow_grid
grid_names=['x' 'y' 'z' 't'];
ref_size = [length(reference.flow_grid.(grid_names(1))) ...
    length(reference.flow_grid.(grid_names(2))) ...
    length(reference.flow_grid.(grid_names(3)))];
for it_grid=1:length(grid_names)-1
    if ~issame(ref_size,size(Velocities.(grid_names(it_grid))))
        error(['The data dimension "' grid_names(it_grid) '" does not match the reference flow_grid.' ]);
    end
end


%-- Convert data to polar (magnitude and angle)
Velocities.magnitude=sqrt(Velocities.x.^2+Velocities.y.^2+Velocities.z.^2);
Velocities.theta=atan2d(-1*Velocities.x,Velocities.z)+180;
ref_Velocities.magnitude=sqrt(reference.Velocities.x.^2+reference.Velocities.y.^2+reference.Velocities.z.^2);
ref_Velocities.theta=atan2d(-1*reference.Velocities.x,reference.Velocities.z)+180;

%-- Apply mask to data
mask_idx=find(ref_Velocities.magnitude == 0);
Velocities.magnitude(mask_idx)=NaN;
Velocities.phi(mask_idx)=NaN;
Velocities.theta(mask_idx)=NaN;

%-- Estimate the error
v_err=abs(Velocities.magnitude-ref_Velocities.magnitude);
angle_err=exp(1i*deg2rad(Velocities.theta))+exp(-1i*deg2rad(ref_Velocities.theta));  % Map to unit circle and substract phases
angle_err=abs(real(angle_err))+1i*(imag(angle_err));  % Make it two quadrants only
angle_err= abs(rad2deg(2*angle(angle_err)));

%-- Make the statitistics relative to the peak value for percentage
v_bias=nanmean(v_err(:))/max(ref_Velocities.magnitude(:))*100;
v_std=nanstd(v_err(:))/max(ref_Velocities.magnitude(:))*100;
ang_bias=nanmean(angle_err(:));        % Relative to 360 degrees
ang_std=nanstd(angle_err(:));

%-- Display scores
fprintf('--------------------------------------------------\n');
fprintf('Velocity magnitude metrics [pct respect %f]: \n',max(ref_Velocities.magnitude(:)));
fprintf('Error: %f, std: %f \n', v_bias, v_std);
fprintf('Velocity angle metrics [deg]: \n');
fprintf('Error: %f, std: %f \n', ang_bias, ang_std);
fprintf('--------------------------------------------------\n');




%% Display figures
set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultfigureposition',[0, 0, 1200, 400]);   % Set aspect ratio
path_output_log=[path_output_results filesep filename_results];

v_limits_plot=[5e-4 max(ref_Velocities.magnitude(:))];
% First figure. The reference CFD vs the estimated
figure();
subplot(1,2,1)
plot_vfi(reference.flow_grid.x*1000,reference.flow_grid.z*1000,squeeze(ref_Velocities.theta)',squeeze(ref_Velocities.magnitude)', ...
    'v_limit',v_limits_plot)
axis image
xlim([-10 10])
xlabel('[mm]')
ylabel('[mm]')
title ('CFD')
% Plot the estimated
subplot(1,2,2)
plot_vfi(reference.flow_grid.x*1000,reference.flow_grid.z*1000,squeeze(Velocities.theta)',squeeze(Velocities.magnitude)', ...
    'v_limit',v_limits_plot)
axis image
xlim([-10 10])
xlabel('[mm]')
ylabel('[mm]')
title ('Ultrasound')
fig = gcf;
fig.PaperPositionMode = 'auto';
if flag_save
    disp(['VFI image saved in "',[path_output_log '_VFI'], '.png','"'])
    print([path_output_log '_VFI'] ,'-dpng','-r0')
end

% Second figure. Plot True values against estimated
a_rad=10;    % circle radius
%%% Sort the datasets
% Magnitude
[true_values,I] = sort(ref_Velocities.magnitude(:));     % Sort the true values
estimated_values =  Velocities.magnitude(:);
% Sort the estimated values by same order
estimated_values = estimated_values(I);
% Sort the angles the same manner
visualization_angles = Velocities.theta(:);
visualization_angles = visualization_angles(I);

% Angles
[true_angles,I] = sort(ref_Velocities.theta(:)); % Sort the true values
estimated_angles =  Velocities.theta(:);
% Sort the estimated values by same order
estimated_angles = estimated_angles(I);
% Sort magnitudes in the same manner
visualization_values = Velocities.magnitude(:);
visualization_values = visualization_values(I);

figure();
ax1=subplot(1,2,1);
plot(true_values,true_values,'r','LineWidth',3 );
hold on
scatter(true_values,estimated_values,a_rad,visualization_angles,'filled')
hold off
ax1.Color = 'black';
ax1.GridColor=[1 1 1];
xlabel('True velocity [m/s]')
ylabel('Estimated velocity [m/s]')
axis image
xlim([0 true_values(end)])
ylim([0 true_values(end)])
grid on
colormap(ax1,hsv(360))
c1=colorbar;
ylabel(c1,'Angle [deg]')


ax2=subplot(1,2,2);
plot(true_angles,true_angles,'r','LineWidth',3 );
hold on
scatter(true_angles,estimated_angles,a_rad,visualization_values,'filled')
hold off
ax2.Color = 'black';
ax2.GridColor=[1 1 1];
axis image
xlabel('True angle [deg]')
ylabel('Estimated angle [deg]')
xlim([0 360])
ylim([0 360])
grid on
c2=colorbar;
ylabel(c2,'Velocity [m/s]')
colormap(ax2,parula)
fig = gcf;
fig.PaperPositionMode = 'auto';
if flag_save
    disp(['Scatter image saved in "',[path_output_log '_scatter'], '.png','"'])
    print([path_output_log '_scatter'] ,'-dpng','-r0')
end

% Save figure if enabled
if flag_save
    disp(['Results saved in "',path_output_log, '.mat','"'])
    results = struct();
    results.Velocities = Velocities;
%     results.n_emissions = n_emissions;
    save(path_output_log,'-struct','results');
end

disp('Evaluation Done')



