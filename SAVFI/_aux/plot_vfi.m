function [h] = plot_vfi(varargin)
%PLOT_VFI               This function displays the vector flow estimates using a different
%                       color and brightness to denote the angle and magnitude respectively.
%                       The estimates can be overlaid over a second image and arrows can be
%                       plotted on top of the velocity image.
%
%
% USAGE :
% [fig_handle] = 	PLOT_VFI(ang,vel)
%
%    INPUT :
%         ang:    angles in a 2-dimensional MxN matrix
%                 [0-360 deg.]
%         vel:    magnitudes a 2-dimensional MxN matrix
%                 [m/s]
%
% [fig_handle] = 	PLOT_VFI(x,z,ang,vel)
%         x,z:   are vectors specifying the locations of the pixel
%                centers of ang(1,1) and ang(M,N).
%
% [fig_handle] = 	PLOT_VFI(x,z,ang,vel,x_bmode,z_bmode,bmode_img)
%         bmode_img:         log_compressed data for generating the bmode image
%                            [0 -db_max] . Default db_max=60
%         x_bmode,z_bmode:   are vectors specifying the locations of the pixel
%                            centers of bmode_img(1,1) and bmode_img(M,N).
%
% [fig_handle] = 	PLOT_VFI(h,x,z,ang,vel,x_bmode,z_bmode,bmode_img)
%         h:      figure handle of the figure to update.
%
% Additional parameters:
%             'v_limit' - [lim_min lim_max] Maximum and minimum velocity to
%                                           display. Default: all
%             'arrows_enable' - [boolean]  Display arrows display.
%                                          Default: true
%             'mask'       - [MxN] Defining the mask to display
%
% Version 1.0 31/03/15 cavh

%% Defualts   %%
par.v_limit = [];
par.arrows_enable = true;
par.cm_option = 1;
par.mask=[];
par.db_max=30;
par.db_mask_max=12.5;
par.bmode_gray_levels=256;
par.arrows_density=[4;8];  
par.arrow_scale=1.5;
par.start_x_arrow=[];
par.start_z_arrow=[];

%% Remove the optional fields %%
params=fieldnames(par);
rm_vargin=[];
% Check for string and get the next value
for it=1:nargin
    if (ischar(varargin{it}))
        % Check if is a valid option
        if(sum(strcmp(varargin{it},params)))
            % Assing value to parameter
            par.(varargin{it})=varargin{it+1};
            % Add to remove list
            rm_vargin=[rm_vargin;it;it+1];
        else
            warning(['Given option ' varargin{it} ' is not valid']);
            % Add to remove list anyway
             rm_vargin=[rm_vargin;it;it+1];
        end
    end
end
% Remove the options
varargin(rm_vargin)=[];
par.bmode_enable = false;

switch (nargin-length(rm_vargin)),
    case 0,
        error('at least two inputs are needed');
    case 1,
        error('at least two inputs are needed');
    case 2,
        hh=gcf;
        curr_ang=varargin{1};
        curr_vel=varargin{2};
        x_axis=[1 size(curr_ang,1)];
        z_axis=[1 size(curr_ang,2)];
    case 3,
        if ~ishandle(varargin{1})
            error('The first argument needs to be a handle')
        end
        hh = varargin{1};
        curr_ang=varargin{1};
        curr_vel=varargin{2};
        x_axis=[1 size(curr_ang,1)];
        z_axis=[1 size(curr_ang,2)];
    case 4,
        hh=gcf;
        x_axis=varargin{1};
        z_axis=varargin{2};
        curr_ang=varargin{3};
        curr_vel=varargin{4};
    case 5,
        if ~ishandle(varargin{1})
            error('The first argument needs to be a handle')
        end
        hh = varargin{1};
        x_axis=varargin{2};
        z_axis=varargin{3};
        curr_ang=varargin{4};
        curr_vel=varargin{5};
    case 6,
        error('The number of input arguments is incorrect');
    case 7,
        hh=gcf;
        x_axis=varargin{1};
        z_axis=varargin{2};
        curr_ang=varargin{3};
        curr_vel=varargin{4};
        x_axis_bmode=varargin{5};
        z_axis_bmode=varargin{6};
        bmode_img=varargin{7};
        par.bmode_enable = true;
    case 8,
        hh = varargin{1};
        x_axis=varargin{2};
        z_axis=varargin{3};
        curr_ang=varargin{4};
        curr_vel=varargin{5};
        x_axis_bmode=varargin{6};
        z_axis_bmode=varargin{7};
        bmode_img=varargin{8};
        par.bmode_enable = true;
    otherwise,
        error('The number of input arguments is incorrect');
end
if (par.bmode_enable)
% Set mask from B-mode
  % Obtain the bmode area that includes the vfi
    x_temp=linspace(x_axis(1),x_axis(end),size(curr_ang,2));
    z_temp=linspace(z_axis(1),z_axis(end),size(curr_ang,1));
    x_temp_bmode=linspace(x_axis_bmode(1),x_axis_bmode(end),size(bmode_img,2));
    z_temp_bmode=linspace(z_axis_bmode(1),z_axis_bmode(end),size(bmode_img,1));
    [x_interp,z_interp] = meshgrid(x_temp,z_temp);
    [x_temp_bmode,z_temp_bmode] = meshgrid(x_temp_bmode,z_temp_bmode);
    bmode_flow=interp2(x_temp_bmode,z_temp_bmode,bmode_img,x_interp,z_interp);
    bmode_flow = medfilt2(bmode_flow, [9 9]);
    mask_bmode=ones(size(bmode_flow));
    mask_bmode(bmode_flow>-par.db_mask_max)=0;
    mask_bmode=mask_bmode;
end

if (isempty(par.mask) && par.bmode_enable)
  

    mask=mask_bmode;
    
elseif (isempty(par.mask))
    mask=ones(size(curr_ang));
elseif (~isempty(par.mask) && par.bmode_enable)
    mask=par.mask.*mask_bmode;
else
    mask=par.mask;
end

if (isempty(par.v_limit))
    par.v_limit=[min(curr_vel(:)) max(curr_vel(:))];
end


% Generate map
map=create_colormap(par.cm_option);

% Limit to 360 and 0
curr_ang(curr_ang>360)=360;
curr_ang(curr_ang<0)=0;
% Magnitude should be only positive
curr_vel(curr_vel<0)= 0;
% Make RGB image of angles
curr_image=round(curr_ang)+1;
rgb_img = ind2rgb(curr_image.*mask,map);
hsv_img = rgb2hsv(rgb_img);
% Change the saturation and value given the velocity (linear scale)
vel_scale=(curr_vel-par.v_limit(1))./(par.v_limit(end)-par.v_limit(1));
vel_scale(vel_scale>1)=1;
vel_scale(vel_scale<0)=0;
hsv_img(:,:,3)=vel_scale.*mask; % Quadratic mapping
% Double check de limits
hsv_img(hsv_img>1)=1;
hsv_img(hsv_img<0)=0;
% Transform back to RGB
rgb_img=hsv2rgb(hsv_img);
rgb_img(rgb_img>1)=1;
rgb_img(rgb_img<0)=0;

figure(hh);

% Update Bmode
if par.bmode_enable
    bmode_img_log=bmode_img+par.db_max;
    bmode_img_log(bmode_img_log<0)=0;
    bmode_img_log=round(bmode_img_log/par.db_max*(par.bmode_gray_levels-1));
    bmode_rgb = ind2rgb(bmode_img_log,gray(par.bmode_gray_levels));
    imagesc(x_axis_bmode,z_axis_bmode,bmode_rgb)
    hold on
else
    % Use a black background
    imagesc([x_axis(1) x_axis(end)],[z_axis(1) z_axis(end)],zeros(size(rgb_img)));
    hold on
end

% Display flow
h=imagesc([x_axis(1) x_axis(end)],[z_axis(1) z_axis(end)],rgb_img);
% Set transparency
curr_mask=mask;
curr_mask(curr_vel<par.v_limit(1))=0;
set(h, 'AlphaData', curr_mask);

% Generate arrows 
if par.arrows_enable
    hold on;
    % Mask angles
    quiv_vel=curr_vel(1:par.arrows_density(1):end,1:par.arrows_density(end):end);
    quiv_angles=curr_ang(1:par.arrows_density(1):end,1:par.arrows_density(end):end);
    % Remove any velocity bigger than the limit
    quiv_vel(quiv_vel>par.v_limit(end))=par.v_limit(end);
    % Scale the values
    quiv_vel= quiv_vel/par.v_limit(end);
    % Make the arrows grid
    x_temp=linspace(x_axis(1),x_axis(end),size(curr_ang,2));
    z_temp=linspace(z_axis(1),z_axis(end),size(curr_ang,1));
    [x,z] = meshgrid(x_temp(1:par.arrows_density(end):end),z_temp(1:par.arrows_density(1):end));
    %Insert a reference arrow on defined point
    quiv_vel_ref=par.v_limit(end);
    quiv_angles_ref=270;
    % Outside mask
    temp_mask=curr_mask(1:par.arrows_density(1):end,1:par.arrows_density(end):end);
    idx_mask=find(temp_mask==1);
    % Get components
    quiv_vx=quiv_vel(idx_mask).*sind(quiv_angles(idx_mask));
    quiv_vz=(-1)*quiv_vel(idx_mask).*cosd(quiv_angles(idx_mask));
    % Add reference arrow
    if (~isempty(par.start_x_arrow)) && (~isempty(par.start_z_arrow))
        quiv_vx=[quiv_vel_ref.*(-1)*sind(quiv_angles_ref); quiv_vx];
        quiv_vz=[(-1)*quiv_vel_ref.*cosd(quiv_angles_ref); quiv_vz];
    end
    x=[par.start_x_arrow; x(idx_mask)];
    z=[par.start_z_arrow; z(idx_mask)];
    quiver(x,z ,...
            quiv_vx*par.arrow_scale, ...  
            quiv_vz*par.arrow_scale, 'w', 'LineWidth', 1,'Autoscale','off')
 
   
end

hold off
end

function map=create_colormap(option)

switch (option),
    % HSV colormap with black 
    case 1,
        map=hsv(360);
        map=[0 0 0; map];
        
    otherwise,
        error('The number of input arguments is incorrect');
end

end


