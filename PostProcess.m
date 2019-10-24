classdef PostProcess < handle
    %PostProcess Container for handling post processing / regularisation of
    %PIV data
    
    properties (Access = public)
        %data
        dTheta
        dR
        Theta
        R
        Wcorr_polar
        vectype_polar
        U
        V
        X
        Z
        x
        z
        Wcorr
        vectype
        removed_outliers
        vdims
        dt    % time between frames
        vdt   % time between vector frames (after correlation averaging)
        t     % time axis
        %copies of scan data - to build matching vector axes
        origin
        ROI
        theta_res
        r_res
        z_res
        x_res
        %options
        do_outlier_removal
        do_universal_outlier
        do_outlier_filling
        do_smoothing
        do_gaussian_spatial_smoothing
        do_temporal_moving_average
        do_smoothn
        %settings
        epsilon
        thresh
        b
        minNXcorrValue
        spatial_smooth_params
        temporal_mav
        is_spatial_disp
        is_spatial_vel
        is_polar
    end
    
    %% constructor
    methods (Access = public)
        function self=PostProcess(piv, scan)
            if strcmp(scan.coord, 'polar')
                self.dTheta = real(piv.U);
                self.dR = real(piv.V);
                self.Theta = real(piv.X);
                self.R = real(piv.Z);
                self.Wcorr_polar = piv.corrMap;
                self.vectype_polar = piv.vflags;
                self.is_polar = true;
                self.theta_res = scan.dtheta;
                self.r_res = scan.dr;
                self.vdims = size(self.dTheta);
            elseif strcmp(scan.coord,'cartesian')
                self.U = real(piv.U);
                self.V = real(piv.V);
                self.X = real(piv.X);
                self.Z = real(piv.Z);
                self.Wcorr = piv.corrMap;
                self.vectype = piv.vflags;
                self.is_polar = false;
                self.vdims = size(self.U);
                
               self.x = scan.x;
               self.z = scan.z;
            end
            self.x_res = scan.dx;
            self.z_res = scan.dz;
            self.origin = scan.origin;
            self.ROI = scan.ROI; 
            self.removed_outliers = zeros(self.vdims(3),1);
            self.dt = scan.dt; 
            self.vdt = Opts.n_ave * self.dt;
            self.t = single(scan.t(1:piv.piv_stop));
            % options
            self.do_outlier_removal = Opts.remove_outliers;
            self.do_universal_outlier = Opts.use_universal_outlier;
            self.do_outlier_filling = Opts.fill_outliers;
            self.do_smoothn = Opts.useSmoothN_robust;
            self.do_smoothing = Opts.smoothing;
            self.do_gaussian_spatial_smoothing = true;
            self.do_temporal_moving_average = true;
            % universal outlier detector settings
            self.epsilon = 0.2;
            self.thresh = 3;
            self.b = 1;
            % other outlier settings
            self.minNXcorrValue = Opts.mincorr; %minimum xcorr value - set to zero if less than
            % smoothing settings
            self.spatial_smooth_params = Opts.sm; % gaussian parameters for spatial smoothing
            self.temporal_mav = Opts.tm; % temporal moving average length
            % initialisation of grid data
            self.is_spatial_disp = false;  % displacement data initialized in pixel units
            self.is_spatial_vel = false; % displacement data initialized in pixels/frame (not m/s)
        end
        
        
        
        %% methods
        function postprocess(self)
            %%% postprocessing of piv results
            if self.do_outlier_removal
                self.remove_outliers();
            end
            if self.do_outlier_filling
                self.fill_outliers();
            end
            if self.do_smoothing
                self.smooth_results();
            end
        end
        
        function set_PIV_results(self,dHorz, dVert)
            %%% sets the applicable PIV matrices depending on if data is
            %%% polar or cartesian
            %%% Horz, Vert are vector positions in polar/cart
            %%% dHorz, dVert are vector component lengths in polar/cart
            if self.is_polar
                self.dTheta = dHorz;
                self.dR = dVert;
            else
                self.U = dHorz;
                self.V = dVert;
            end
        end
        
        function [dHorz, dVert, Corr, VMask] = get_PIV_results(self)
            %%% gets the applicable PIV matrices depending on if data is
            %%% polar or cartesian
            %%% dHorz, dVert are vector component lengths in polar/cart
            %%% Corr matrix stack of maximum normalized cross correlation values per
            %%% kernel
            %%% VMask is matrix stack indicating if vector is masked
            if self.is_polar
                dHorz = self.dTheta;
                dVert = self.dR;
                Corr = self.Wcorr_polar;
                VMask = self.vectype_polar;
            else
                dHorz = self.U;
                dVert = self.V;
                Corr = self.Wcorr;
                VMask = self.vectype;
            end
        end
        
        function remove_outliers(self)
            %%% removes outliers from velocity results
            % first get applicable velocity matrices (polar/cart)
            [dHorz, dVert, Corr, VMask] = self.get_PIV_results();
            % remove lower threshold of correlation values
            dHorz(Corr <self.minNXcorrValue) = nan;
            dVert(Corr <self.minNXcorrValue) = nan;
            % remove flagged values
            dHorz(VMask == -1) = nan;
            dVert(VMask == -1) = nan;
            if self.do_universal_outlier
                % local median outlier (universal outlier - Westerweel)
                for i = 1:size(dHorz,3)
                    [tempu,tempv, self.removed_outliers(i)] = median_outlier(dHorz(:,:,i),...
                        dVert(:,:,i),self.epsilon,self.thresh,self.b);
                    
                    dHorz(:,:,i) = tempu;
                    dVert(:,:,i) = tempv;
                end
            end
            % now set the correct velocity matrices
            self.set_PIV_results(dHorz, dVert);
        end
        
        function fill_outliers(self)
            %%% fills all nan values by interpolating by surrounding data
            % first get applicable velocity matrices (polar/cart)
            [dHorz, dVert, Corr, ~] = self.get_PIV_results();
            % first set kernels with low nxcorr values to 0. This assumes
            % that only noise was present after tissue/clutter suppression
            dHorz(Corr <self.minNXcorrValue) = 0;
            dVert(Corr <self.minNXcorrValue) = 0;
            % now fill all nan values with surrounding data
            for i = 1:size(dHorz,3)
                dHorz(:,:,i)=single(inpaint_nans(double(dHorz(:,:,i)),4));
                dVert(:,:,i)=single(inpaint_nans(double(dVert(:,:,i)),4));
            end
            % now set the correct velocity matrices
            self.set_PIV_results(dHorz, dVert);
        end
        
        function smooth_results(self)
            %%% smooths PIV data
            % first get applicable velocity matrices (polar/cart)
            [dHorz, dVert, ~, VMask] = self.get_PIV_results();
            % first set masked values to 0 so they dont destroy other data
            dHorz(VMask <= 0) = 0;
            dVert(VMask <= 0) = 0;
            if self.do_smoothn
                self.smoothnRobust();
            elseif self.do_gaussian_spatial_smoothing
                self.gauss_spatial_smoothing();
            end
            if self.do_temporal_moving_average
                if self.vdims(3) > self.temporal_mav
                    self.temp_move_ave();
                else
                    warning('Chosen value for temporal_mav is larger than the number of vector frames, ignoring temporal moving average');
                end
            end
            % set masked values back to nan
            dHorz(VMask == 0) = nan;
            dVert(VMask == 0) = nan;
            % now set the correct velocity matrices
            self.set_PIV_results(dHorz, dVert);
        end
        
        function smoothnRobust(self)
           %%% calls Damien Garcia's smoothn function using the robust smoothing mode
           [dHorz, dVert, ~, ~] = self.get_PIV_results();
           vel = {dHorz,dVert};
           vel = smoothn(vel,'robust');
           dHorz = vel{1};
           dVert = vel{2};
           % now set the correct velocity matrices
           self.set_PIV_results(dHorz, dVert);
            
        end
     
        function gauss_spatial_smoothing(self)
            %%% gaussian kernel smoothing in 2D (spatial only)
            % first get applicable velocity matrices (polar/cart)
            [dHorz, dVert, ~, ~] = self.get_PIV_results();
            sm = self.spatial_smooth_params;
            for i = 1:size(self.dTheta,3)
                dHorz(:,:,i)=imgaussfilt(dHorz(:,:,i),[sm(1), sm(2)],...
                    'FilterSize', [sm(3), sm(4)], 'FilterDomain','spatial');
                dVert(:,:,i)=imgaussfilt(dVert(:,:,i),[sm(1), sm(2)],...
                    'FilterSize', [sm(3), sm(4)], 'FilterDomain','spatial');
            end
            % now set the correct velocity matrices
            self.set_PIV_results(dHorz, dVert);
        end
        
        function temp_move_ave(self)
            %%% temporal moving average
            % first get applicable velocity matrices (polar/cart)
            [dHorz, dVert, ~, ~] = self.get_PIV_results();
            dHorz = movingmean(dHorz,self.temporal_mav,3,1);
            dVert = movingmean(dVert,self.temporal_mav,3,1);
            % now set the correct velocity matrices
            self.set_PIV_results(dHorz, dVert);
        end
        
        function convert_to_spatial(self)
            %%% converts vector maps from pixel displacements to
            %%% displacement in mm (rad for polar azimuth) (does not convert to velocities yet)
            if ~self.is_spatial_disp  % only convert to spatial if not already spatial
                if self.is_polar
                    % vector bases in polar domain (need to add ROI offset and
                    % half of kernel size offset to make vector base in
                    % center of kernel)
                    self.Theta = (self.Theta + self.ROI(1)).*self.theta_res...
                        - self.origin(2); % [rad]
                    self.R = (self.R + self.ROI(2)).*self.r_res; % [mm]
                    % vector lengths in polar domain
                    self.dTheta = self.dTheta.*self.theta_res; % [mm]
                    self.dR = self.dR.*self.r_res; % [mm]
                else
                    % vector bases in cartesian domain (need to add ROI offset and
                    % half of kernel size offset to make vector base in
                    % center of kernel)
                    self.X = (self.X + self.ROI(1)).*self.x_res + self.origin(2); % [mm]
                    self.Z = (self.Z + self.ROI(2)).*self.z_res + self.origin(1); % [mm]
                    % vector lengths in cartesian domain
                    self.U = self.U.*self.x_res; % [mm]
                    self.V = self.V.*self.z_res; % [mm]
                end
                self.is_spatial_disp = true;  % set spatial flag so we dont convert to spatial twice
            else
                warning('Data is already in spatial units, ignoring');
            end
        end
        
        function scan_convert_vectors(self)
            %%% converts vectors from polar coordinate system to cartesian
            %%% coordinate system
            if self.is_polar
                if ~self.is_spatial_disp
                    self.convert_to_spatial()
                end
                %tips of vectors
                tip_theta = self.Theta + self.dTheta;
                tip_r = self.R + self.dR;
                %now scan convert bases
                self.X = self.R.*sin(self.Theta); %mm
                self.Z = self.R.*cos(self.Theta); %mm
                %now scan convert tips
                x_tip = tip_r.*sin(tip_theta); %mm
                z_tip = tip_r.*cos(tip_theta); %mm
                %subtract tip from base to get scan converted vectors
                self.U = (x_tip - self.X); % mm
                self.V = (z_tip - self.Z); % mm
                self.Z = self.Z + self.origin(1); % offset by virtual focus
            else
                warning('Tried to scan convert non-polar data, ignoring')
            end
        end
        
        function scan_convert_grid(self, grid_res)
            %%% converts polar grid to cartesian grid for vector data
            %%% grid res [dx, dz] - grid spacing in z and x direction
            %define cartesian grid
            assert(self.is_polar, 'data is not polar, cannot convert from polar to cartesian');
            assert(self.is_spatial_disp,'First convert vectors to spatial coordinate system using convert_to_spatial')
            if nargin < 2  % defaults
                %res in z/x directions (x defined at mean image depth)
                mid_index = fix(size(self.X)./2);
                meandX = mean(diff(self.X,1,2),1);
                meandZ = mean(diff(self.Z,1,1),2);
                grid_res = [meandZ(mid_index(1)), meandX(mid_index(2))]; % dz, dx
            end
            xlims = (max(self.R(:))).*sin([min(self.Theta(:)), max(self.Theta(:))]); %min, max in [mm]
            minz = min(min(self.R(:)).*cos([min(self.Theta(:)), max(self.Theta(:))]));
            maxz = max(max(self.R(:)).*cos([min(self.Theta(:)), max(self.Theta(:))]));
            zlims = [minz, maxz]; %min, max in [mm]
            % define new grid
            [self.Z, self.X] = ndgrid(zlims(1):grid_res(1):zlims(2), xlims(1):grid_res(2):xlims(2));
            % convert vector mask
            self.vectype = uint8(scan_convert(single(self.vectype_polar),self.R, self.Theta, self.Z, self.X));
            % convert maximum nxcorr map
            self.Wcorr = scan_convert(self.Wcorr_polar,self.R, self.Theta, self.Z, self.X);
            % now scan convert vectors
            if self.is_spatial_vel % if data is velocity then first convert to disp
                self.convert_vel_to_disp();
            end
            self.U = scan_convert(self.U, self.R, self.Theta, self.Z, self.X);
            self.V = scan_convert(self.V, self.R, self.Theta, self.Z, self.X);
            self.Z = self.Z + self.origin(1); % offset by virtual focus
        end
        
        function convert_disp_to_vel(self)
            assert(~isempty(self.dt), 'first set vector temporal resolution using set_dt');
            if ~self.is_spatial_vel  % dont convert to velocity twice
                if ~self.is_spatial_disp % if displacements are still in pixel units, first convert to spatial
                    self.convert_to_spatial();
                end
                self.U = self.U./(self.dt*Opts.frameSkip*1e3);  % m/s
                self.V = self.V./(self.dt*Opts.frameSkip*1e3);  % m/s
                self.is_spatial_vel = true;
            else
                warning('Tried converting velocity data to velocity units, ignoring');
            end
            
        end
        
        function convert_vel_to_disp(self)
            assert(~isempty(self.dt), 'first set vector temporal resolution using set_dt');
            if self.is_spatial_vel
                self.U = self.U.*(self.dt*Opts.frameSkip*1e3);  % mm
                self.V = self.V.*(self.dt*Opts.frameSkip*1e3);  % mm
                self.is_spatial_vel = false;
            else
                warning('Tried converting displacement data to displacment units, ignoring')
            end
        end
        
        function set_dt(self, dt)
            %%%  sets the temporal resolution for converting between
            %%%  displacement and velocity units
            self.dt = dt;
        end
        
        function save_piv_data(self, filepath, scan)
            %%% saves PIV data and bmode images (cartesian only)
            %%% filepath - savepath without extensions
            %%% scan - instance of ImageProperties class
            [BmodeShape, BmodeType, bmode_path] = scan.save_bmode_bin(filepath);
            mask = scan.imMask_cart;
            roi = scan.ROI;
            xfilt = self.X;  %mm
            zfilt = self.Z;  %mm
            ufilt = self.U;  %m/s
            vfilt = self.V;  %m/s
            vmaskfilt = self.vectype;
            tRes = self.vdt; % for display we need time between vector frames
            origin = [scan.x(1), scan.x(end), scan.z(1), scan.z(end)]; % image extents
            Wcorrfilt = self.Wcorr;
            unit = 'm/s';
            n_ave = Opts.n_ave;
            xRes = scan.dx;
            zRes = scan.dz;
            piv_start = scan.start_frame;
            piv_stop = scan.end_frame;
            mPoly =  scan.mPoly;
            savepath = [filepath '_piv.mat'];
            save(savepath,'roi','mask','ufilt', 'vfilt','xfilt','zfilt','xRes','zRes',...
                'tRes','origin','piv_start','piv_stop','vmaskfilt','Wcorrfilt',...
                'BmodeShape','BmodeType','bmode_path','mPoly','n_ave','unit','-v6');
        end
          
        function save_SAVFI_results(self, savepath, n_emissions)
            %%% saves SAVFI data for use with the SAVFI evaluator
            %%% filepath = path to original beamformed SAVFI file
            %%% savepath = path to save destination
            %%% n_emissions largest ensemble size used for averaging/estimation
            assert(strcmp(Opts.datasource, 'SAVFI'), 'Opts.datasource not set to SAVFI')
            assert(self.is_spatial_vel, 'First convert vector data to m/s');
            
            flow_grid = get_SAVFI_grid(savepath);
            zax = self.Z(:,1)*1e-3;  % convert to m
            xax = self.X(1,:)*1e-3;  % convert to m
            tax = self.t;
            tmaxidx = self.find_closest_idx(flow_grid.t, nanmax(tax(:)));
            flow_grid.t = flow_grid.t(1:tmaxidx);
            % Reference grid
            Velocities = struct();  % Velocity information
            Velocities.x = zeros(length(flow_grid.z),length(flow_grid.x), tmaxidx);
            Velocities.y = zeros(length(flow_grid.z),length(flow_grid.x), tmaxidx);
            Velocities.z = zeros(length(flow_grid.z),length(flow_grid.x), tmaxidx);
            % interpolate onto reference grid
            fint = griddedInterpolant({zax, xax', tax'},self.U,'linear');
            Velocities.x = fint({flow_grid.z, flow_grid.x, flow_grid.t});
            fint.Values = self.V;
            Velocities.z = fint({flow_grid.z, flow_grid.x, flow_grid.t});
            for comp = ['x', 'y', 'z']
                Velocities.(comp) = permute(Velocities.(comp), [2, 4, 1, 3]);
            end
            % now set nan/masked sections to 0 as eval script does not do
            % so good with nans
            mag = sqrt(Velocities.x.^2 + Velocities.y.^2 + Velocities.z.^2);
            Velocities.x(isnan(mag)) = 0;
            Velocities.y(isnan(mag)) = 0;
            Velocities.z(isnan(mag)) = 0;
            % save data to file so it can be loaded via the evaluation
            % script
            savepath = [savepath '_SAVFI.mat'];
            save(savepath,'Velocities','n_emissions');          
        end
        
        function get_SAVFI_results(self, savepath, n_emissions)
            %%% gets SAVFI results data for use with the SAVFI evaluator
            %%% does not save Velocities structure to file, instead calls
            %%% the applicable evaluation function and saves the results of
            %%% that directly.
            %%% savepath = path to save destination
            %%% n_emissions largest ensemble size used for averaging/estimation
            assert(strcmp(Opts.datasource, 'SAVFI'), 'Opts.datasource not set to SAVFI')
            assert(self.is_spatial_vel, 'First convert vector data to m/s');
            addpath('SAVFI');
            flow_grid = get_SAVFI_grid(savepath);
            zax = self.Z(:,1)*1e-3;  % convert to m
            xax = self.X(1,:)*1e-3;  % convert to m
            tax = self.t;
            tmaxidx = self.find_closest_idx(flow_grid.t, nanmax(tax(:)));
            flow_grid.t = flow_grid.t(1:tmaxidx);
            % Reference grid
            Velocities = struct();  % Velocity information
            Velocities.y = zeros(length(flow_grid.z),length(flow_grid.x), tmaxidx); % 0 out of plane motion
            % interpolate onto reference grid
            fint = griddedInterpolant({zax, xax', tax'},self.U,'linear');
            Velocities.x = fint({flow_grid.z, flow_grid.x, flow_grid.t});
            fint.Values = self.V;
            Velocities.z = fint({flow_grid.z, flow_grid.x, flow_grid.t});
            for comp = ['x', 'y', 'z']
                Velocities.(comp) = permute(Velocities.(comp), [2, 4, 1, 3]);
            end
            % now set nan/masked sections to 0 as eval script does not do
            % so good with nans
            mag = sqrt(Velocities.x.^2 + Velocities.y.^2 + Velocities.z.^2);
            Velocities.x(isnan(mag)) = 0;
            Velocities.y(isnan(mag)) = 0;
            Velocities.z(isnan(mag)) = 0;
            if strfind(savepath, 'straight_vessel') ~= -1
                evaluation_straight_vessel(Velocities, n_emissions, savepath)
            elseif strfind(savepath, 'carotid_bifurcation') ~= -1
                evaluation_cfd_model(Velocities, savepath)
            elseif strfind(savepath, 'spinning_disk') ~= -1
                evaluation_spinning_disk(Velocities, n_emissions, savepath)
            else
                error('Could not determine the phantom type from the path')
            end
                
                
        end
        
        function [I] = find_closest_idx(~,arr, val)
            %%% finds the index of the closest value in arr to val
            [~,  I] = nanmin(abs(arr-val));
        end
            
        function play_velmag(self, lims)
            %%% plays movie of velocity magnitudes for quick debugging
            %%% lims = [min max] velocity magnitude for display
            figure()
            hold on
            mag = sqrt(self.U.^2 + self.V.^2);
            mag = flipud(mag);
            %hax = [nanmin(self.X(:)), nanmax(self.X(:))];
            hax = [self.x(1), self.x(end)];
            vax = [self.z(1), self.z(end)];
            %vax = [nanmin(self.Z(:)), nanmax(self.Z(:))];
            if nargin < 2
                lims = [prctile(mag(:), 1) prctile(mag(:), 99.9)];
            end
            fig = imagesc(hax,vax, mag(:,:,1));%,lims);
            colormap hot
            axis image
            c = colorbar();
            if self.is_spatial_vel
                c.Label.String = 'Velocity Magnitude [m/s]';
            else
                c.Label.String = 'Displacement Magnitude [mm]';
            end
            xlabel('Width [mm]')
            ylabel('Depth [mm]')
            for i = 1:size(mag,3)
                fig.CData = mag(:,:,i);
                pause(0.01)
            end
        end
        
    end % methods
end %class

