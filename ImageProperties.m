classdef ImageProperties < handle
    % ImageProperties axes / scan / resolution properties and functions
    % of images used for PIV
    properties (Access = public)
        % begin necessary info
        filepath      % path to where data is stored
        frames        % ultrasound image sequence - can be cart or polar (r/z, theta/x, angles, time)
        % if cartesian
            dx            % x resolution (cartesian - width) [mm]
            dz            % z resolution (cartesian - depth) [mm]
            x             % x axis [mm] (cartesian - width)
            z             % z axis [mm] (cartesian - depth)
        % if polar
            dtheta        % theta pixel resolution in radians
            dr            % radial pixel resolution in mm
            r             % depth axis [mm]
            theta         % azimuthal axis [rad]
        origin        % origin of frame coordinate system (vert, horizontal)
        %either
            dt            % temporal resolution (s) between like tx/rx acqs.
        % or
            FR            % frame rate between like tx/rx acqs.
        %end necessary information
        bmodes_polar  % frames for showing vectors on in polar domain
        bmodes_cart   % frames for showing vectors on in cartesian domain
        frame_size    % size of frame (pixels)
        num_ang       % number of angles
        num_frames    % number of time points
        num_frames_bmode % number of frames for bmode (may be different due to ensemble averaging)
        start_frame   % frame to start PIV analysis
        end_frame     % frame to end PIV analysis
        fs            % VSX operating frequency
        cf            % received pulse center frequency
        c             % speed of sound [m/s]
        coord         % coordinate system data is defined in ('polar', 'cartesian') (default 'cartesian');
        z_0           % z axis [mm] (cartesian - depth - referenced to virtual focus) [mm]
        t             % time axis [s]
        scan_mode     % specify whether to scan convert whole image or part of ('full', 'cropped')
        %% mask properties
        ROI           % region of interest (pixels)
        ROI0          % ROI unmodified by padding
        imMask_polar        % boolean mask indicating masked regions with true
        imMask_cart   % boolean mask scan converted
        mPoly         % hand drawn polygon
        pad_inside    % padding on left top side of image (to fit kernel)
        pad_outside   % padding on right bottom side of image (to fit kernel)
        %% visualisation
        log_compression % quick function for taking envelope and log compressing
        
    end
    %% constructor
    methods (Access = public)
        function self=ImageProperties()
            self.log_compression = @(x) 20*log10(abs(x/max(x(:))));  % handy function for log compression
        end
    end
    %% methods
    methods
        function update(self)
            %%% update calculated parameters after the necessary parameters
            %%% have been manually filled
            assert(~isempty(self.frames), 'need frames in shape (r/z, theta/x, angles, time)')
            if isempty(self.coord)  % first establish coordinate system used
                self.coord = 'cartesian';
            end
            assert(ndims(self.frames) == 4, 'frames should have shape (r/z, theta/x, angles, time)') 
            fsize = size(self.frames);  % get data dimensions        
            self.frame_size = fsize(1:2);
            self.num_ang = fsize(3);
            self.num_frames = fsize(4);
            
            if strcmp(self.coord, 'cartesian')
                assert(~(isempty(self.dx) || isempty(self.dz) || isempty(self.origin))...
                    && (isempty(self.x) || isempty(self.z)), ...
                    'Need either dx/dz and frames or x and z axes')
                if isempty(self.x) || isempty(self.z)  % supplied origin, dx and dz
                    self.x = self.origin(2) + [0:self.frame_size(2)].*self.dx;
                    self.z = self.origin(1) + [0:self.frame_size(1)].*self.dz;
                else  % supplied x and z axes
                    self.origin = [self.z(1), self.x(1)];
                    self.dx = mean(diff(self.x,1));
                    self.dz = mean(diff(self.z,1));
                end
            else  % polar coordinate system - TODO: untested - origin is probably wrong - wrote quick and dirty
                assert(~(isempty(self.dtheta) || isempty(self.dr) || isempty(self.origin))...
                    && (isempty(self.theta) || isempty(self.r)), ...
                    'Need either dtheta/dr and frames or theta and r axes')
                if isempty(self.theta) || isempty(self.r)  % supplied origin, dtheta and dr
                    self.theta = self.origin(2) + [0:self.frame_size(2)].*self.dtheta;
                    self.r = self.origin(1) + [0:self.frame_size(1)].*self.dr;
                else  % supplied theta and r axes
                    self.origin = [self.r(1), self.theta(1)];
                    self.dtheta = mean(diff(self.theta,1));
                    self.dr = mean(diff(self.r,1));
                end
            end
            assert(~(isempty(self.dt) && isempty(self.FR)), 'Need either dt or frame rate')
            if isempty(self.dt)
                self.dt = 1/self.FR; 
            else
                self.FR = 1/self.dt;
            end
            self.t = [0:self.num_frames].*self.dt;
            % get number of frames
            self.get_frame_lims()
            self.trim_frames()  % trim frames we dont want to process
            
        end
        
        function loadFromBinFileWrapper(self, filepath)
            %%% loadFromBinFileWrapper Thin wrapper for loading a binfile
            %%% Load Data from Binfile
            self.filepath = filepath;
            datapath = filepath(1:strfind(filepath, '_1_info.mat')-1);
            datapath = [datapath '.bin'];
            % r/z, theta/x, angles, time
            self.frames = loadBinFile(filepath, datapath);
            load(filepath, 'P', 'PData');  % load P and PData info
            self.getInfoFromVSX(P, PData);  % get the info needed from the vsx structures
            % r/z, theta/x, angles, time
            self.frames = reshape(self.frames, self.frame_size, self.num_ang, self.num_frames);
            self.get_frame_lims()
            self.trim_frames()  % trim frames we dont want to process
            
        end
        
        function loadFromMatfilesWrapper(self, filepath)
            %%% loadFromMatfilesWrapper Loads files saved as V6 matfile series. Need to
            %%% specify the info file
            self.filepath = filepath;
            A = load(filepath);
            P = A.P;
            PData = A.PData;
            self.getInfoFromVSX(P, PData);  % get the info needed from the vsx structures
            batch_size = P.numAcqs*P.na;
            numacqs = batch_size*P.numFrames;
            numframes = numacqs/P.na;
            if (Opts.procFrames==0)||(Opts.procFrames>numframes)
                start_batch = 1;
                end_batch = numacqs/batch_size;
                numbatches = (end_batch - start_batch + 1);
            else
                start_batch = floor(P.na*Opts.startFrame/batch_size)+1;
                end_batch = ceil(P.na*Opts.endFrame/batch_size);
                numbatches = (end_batch - start_batch + 1);
            end
            self.num_frames = numbatches*batch_size/self.num_ang;
            dims = [self.frame_size, self.num_frames*self.num_ang];  %loading in 3 dims only (time*angles for dim 3)
            self.frames = zeros(dims, 'single');
            datapath = filepath(1:strfind(filepath, '_1_info.mat')-1);
            fnameid = cell(numbatches);
            k = 1;
            for i = start_batch:end_batch
                fnameid{k} = [datapath '_' num2str((i-1)*batch_size+1) '.mat'];
                A = load(fnameid{k});
                self.frames(:,:,batch_size*(k-1)+1:k*batch_size) = squeeze(A.IQData);
                k = k+1;
            end
            
            % r/z, theta/x, angles, time
            self.frames = reshape(self.frames, [self.frame_size, self.num_ang, self.num_frames]);
            % get number of frames
            self.get_frame_lims()
            self.trim_frames()  % trim frames we dont want to process
            
        end
        
        function get_frame_lims(self)
            %%% gets the start and stop frame number and number of frames
            %%% to process based on the frames contained in the data and
            %%% the options chosen by the user
            if (Opts.procFrames==0)||(Opts.procFrames>self.num_frames)
                self.start_frame = Opts.startFrame;
                self.end_frame = self.num_frames;
                self.num_frames = self.end_frame - self.start_frame + 1;
            else
                self.num_frames = Opts.procFrames;
                self.start_frame = Opts.startFrame;
                self.end_frame = Opts.endFrame;
            end
        end
        
        function trim_frames(self)
            %%% trims frames to work with the number of averages
            %%% specified in Opts class and start and end frames

            % make sure divisable by num averages
            self.num_frames = self.num_frames - mod(self.num_frames,Opts.n_ave);
            self.end_frame = self.start_frame + self.num_frames - 1;
            
            %trim frames we dont want to process
            % note the frame 1 in frames is already offset by start_frame -
            % only need to trim end
            self.frames = self.frames(:,:,:,1:self.num_frames);
            self.t = ([self.start_frame:Opts.n_ave:self.end_frame]-1).*self.dt; %update time axis [s]
        end
        
        function choose_angles(self, angs)
            %%% convenience function for chosing a subset of available
            %%% angles. Overwrites self.frames data.
            %%% angs - list of indices for angles that are wanting to be
            %%% kept
            self.frames = self.frames(:,:,angs,:);
            self.num_ang = length(angs);
        end
        
        function getInfoFromVSX(self, P, PData)
            %%% getInfoFromVSX converts the beamforming information used by Verasonics
            %%% into an info file compatible with PIV.
            %%%   At present only works with a single PData struct and 2D scans only.
            %%% P - P structure used when beamforming
            %%% PData - VSX PData structure
            self.num_ang = P.na;
            self.FR = 1/(P.PRF*P.na*1e-6);
            try self.fs = P.opFreq; catch; self.fs = P.txFreq; end
            try self.cf = P.centerFreq; catch; self.cf = self.fs; end
            self.c = P.c;
            wvl2mm = (1e-3*self.c/self.fs);  % for computing VSX beamforming grid
            if strcmp(PData.Coord, 'polar') % beamformed onto polar grid
                %no steering accounted for
                self.coord = 'polar';
                self.dtheta = PData.PDelta(1);
                self.dr = PData.PDelta(2)*wvl2mm;
                self.origin = [PData.Origin(3)*wvl2mm, (PData.Size(2)/2)*PData.PDelta(1)];
                try r2 = PData.Region.Shape.r;catch r2 = PData.Region.Shape.r2;end %in wvl
                thetaRange = PData.Region.Shape.angle;
                %         thetaSteer = PData.Region.Shape.steer;  %polar doesnt work with steered regions
                rlims = [1, r2].*wvl2mm; % r_limits of whole image
                thetalims = [self.origin(2) - thetaRange, thetaRange - self.origin(2)]; %min, max in rad VSX limits theta to being symmetric around origin
                self.frame_size = PData.Size(1:2);
                self.r = linspace(rlims(1), rlims(2), self.frame_size(1)); % referenced to apex of sector
                self.theta = linspace(thetalims(1), thetalims(2), self.frame_size(2));
            elseif strcmp(PData.Coord, 'rectangular')
                %no steering accounted for
                self.coord = 'cartesian';
                self.dx = PData.PDelta(1)*wvl2mm; % x pixel resolution in mm
                self.dz = PData.PDelta(2)*wvl2mm;  % z pixel resolution in mm
                self.origin = [PData.Origin(3), PData.Origin(1)].*wvl2mm; % [z, x] position of top left corner pixel in mm
                self.frame_size = PData.Size(1:2);
                xlims = [self.origin(2), self.origin(2)+self.frame_size(2)*self.dx];
                zlims = [self.origin(1), self.origin(1)+self.frame_size(1)*self.dz];
                self.z = linspace(zlims(1), zlims(2), self.frame_size(1)); % z pixel positions [mm]
                self.x = linspace(xlims(1), xlims(2), self.frame_size(2)); % x pixel positions [mm]
            end
            self.dt = 1/self.FR;  % temporal resolution between like angles
            self.t = 0:self.dt:(self.num_frames-1)*self.dt; % time axis [s]
        end
        
        function loadSAVFIdata(self, filepath)
            %%% loads data from the synthetic aperture vector flow imaging
            %%% (SAFVI) challenge dataset. Uses the included beamformer.
            %%% Assumes all beamformed data was save with variable name
            %%% 'Frames' and shape is [z, x, angles, time]
            %%% needs 'x_axis' and 'z_axis' used for beamforming to be
            %%% saved as well as arrays.
            self.filepath = filepath;
            A = load(filepath);
            self.frames = loadComplexBin(A.bin_path, A.SaveShape, A.DataType, A.DataShape);
            dims = size(self.frames);
            self.frame_size = dims(1:2);
            self.num_ang = dims(3);  %number of Sythetic Apertures / Angles
            self.num_frames = dims(4);
            self.x = A.x_axis*1e3; % xaxis - SAVFI is in m we want mm
            self.z = A.z_axis*1e3; % zaxis - SAVFI is in m we want mm
            self.FR = A.para.sys.fprf/(self.num_ang+1);  % frame rate between like 
            % angles, + 1 because of bmode interleved bmode sequence frame
            self.fs = A.para.xdc.f0;  % transmit wave frequency
            self.cf = self.fs;
            self.c = A.para.sys.c;  %speed of sound
            self.coord = 'cartesian';
            self.dx = mean(diff(self.x,1)); % x pixel resolution in mm
            self.dz = mean(diff(self.z,1));  % z pixel resolution in mm
            self.origin = [self.z(1), self.x(1)]; % [z, x] position of top left corner pixel in mm
            self.dt = 1/self.FR;  % temporal resolution between like angles
            self.t = 0:self.dt:(self.num_frames-1)*self.dt; % time axis [s]
            self.get_frame_lims()
            self.trim_frames()  % trim frames we dont want to process
            
        end
        
        function loadSAVFIustbdata(self, filepath)
            %%% loads data from the synthetic aperture vector flow imaging
            %%% (SAFVI) challenge dataset which is saved using USTB as a uff file.
            self.filepath = filepath;
            contents = uff.index(filepath, '/', true);
            % find variable that is of beamformed_data class
            bfidx = strfind(cellfun( @(sas) sas.class, contents, 'uni', false ), {'uff.beamformed_data'});
            bfdata_name = contents{bfidx{1}}.name;
            bfdata = uff.beamformed_data();
            bfdata.read(filepath,['/' bfdata_name]);
            self.frames = bfdata.data;
            dims = size(bfdata.data);
            self.frame_size = dims(1:2);
            self.num_ang = dims(3);  %number of Sythetic Apertures / Angles
            self.num_frames = dims(4);
            self.x = bfdata.scan.x_axis*1e3; % xaxis - SAVFI is in m we want mm
            self.z = bfdata.scan.z_axis*1e3; % zaxis - SAVFI is in m we want mm
            self.FR = bfdata.PRF;  % frame rate between like angles/modes
            self.fs = bfdata.pulse.center_frequency;  % transmit wave frequency
            self.cf = bfdata.pulse.center_frequency;  % assume center freq is same as mod freq
            self.c = bfdata.sequence(1).sound_speed;  %speed of sound
            if contains(class(bfdata.scan), 'linear_scan')
                self.coord = 'cartesian';
            elseif contains(class(bfdata.scan), 'sector_scan')
                self.coord = 'polar';
            else
                error('Could not determine the scan type of the uff file')
            end
            self.dx = mean(diff(self.x,1)); % x pixel resolution in mm
            self.dz = mean(diff(self.z,1));  % z pixel resolution in mm
            self.origin = [self.z(1), self.x(1)]; % [z, x] position of top left corner pixel in mm
            self.dt = 1/self.FR;  % temporal resolution between like angles
            self.t = 0:self.dt:(self.num_frames-1)*self.dt; % time axis [s]
            self.get_frame_lims()
            self.trim_frames()  % trim frames we dont want to process=
        end
        
        
        function printPIVInfo(self)
            %%% printPIVInfo Prints info related to kernel size and tracking
            %%% (currently only works for polar)
            if strcmp(self.coord, 'polar')
                ktheta = Opts.ksize(:,2).*self.dtheta; %rad
                kthetadeg = ktheta*180/pi; %deg
                kr = Opts.ksize(:,1).*self.dr; %mm
                for i = 1:Opts.nIterations
                    fprintf('The kernel size for iteration %d in spatial coords is: %4.2f degrees x % 4.2f mm \n', i, kthetadeg(i), kr(i));
                end
                % maximum trackable velocity according to 1/4 rule
                midr = (self.r(end) - self.r(1))/2; %mm
                midx = midr*sin(ktheta(1)/4); %mm
                maxvel = [midx/self.dt (kr(1)/4)/self.dt]./1e3; %m/s
                fprintf('Maximum theoretical trackable velocity in midpoint of image is %4.1f m/s x %4.1f m/s in azimuth/axial directions \n', maxvel(1), maxvel(2));
                fprintf('Maximum theoretical trackable velocity in midpoint of image, taking into account frameskip, is %4.1f m/s x %4.1f m/s in azimuth/axial directions \n', maxvel(1)/Opts.frameSkip, maxvel(2)/Opts.frameSkip);
                %final grid size
                finalgridxy = [midr*sin(Opts.step(Opts.nIterations,2)*self.dtheta), Opts.step(Opts.nIterations,1)*self.dr];
                finalgridrtheta = [Opts.step(Opts.nIterations,2)*self.dtheta*180/pi, Opts.step(Opts.nIterations,1)*self.dr];
                fprintf('Final grid spacing = %4.2f mm by % 4.2f mm (x/z) \n',finalgridxy(1), finalgridxy(2))
                fprintf('Final grid spacing = %4.2f ? by % 4.2f mm (theta/r) \n',finalgridrtheta(1), finalgridrtheta(2))   
                % smoothing kernel size
                smstd = Opts.sm(1:2).*finalgridxy;
                smext = Opts.sm(3:4).*finalgridxy;
                fprintf('Spatial smoothing kernel: std dev = %4.2f mm x %4.2f mm, extent = %4.2f mm x %4.2f mm \n', smstd, smext)
                fprintf('Temporal smoothing kernel = %4.2f ms\n',Opts.tm*self.dt*1e3)
            else
                kx = Opts.ksize(:,2).*self.dx;  % kernel size along x axis [mm]
                kz = Opts.ksize(:,1).*self.dz;  % kernel size along z axis [mm]
                for i = 1:Opts.nIterations
                    fprintf('The kernel size for iteration %d in spatial coords is: x = %4.2f mm x z = %4.2f mm \n', i, kx(i), kz(i));
                end
                % maximum trackable velocity according to 1/4 rule
                maxvel = [kx(1) kz(1)]/(4*Opts.frameSkip*self.dt*1e3); %m/s
                fprintf('Maximum theoretical trackable velocity, with frameSkip is %4.1f m/s x %4.1f m/s in x/z directions \n', maxvel(1), maxvel(2));
                %final grid size
                finalgridxy = [Opts.step(Opts.nIterations,2)*self.dx, Opts.step(Opts.nIterations,1)*self.dz];
                fprintf('Final grid spacing = %4.2f mm by % 4.2f mm (x/z) \n',finalgridxy(1), finalgridxy(2))
                % smoothing kernel size
                smstd = Opts.sm(1:2).*finalgridxy;
                smext = Opts.sm(3:4).*finalgridxy;
                fprintf('Spatial smoothing kernel: std dev = %4.2f mm x %4.2f mm, extent = %4.2f mm x %4.2f mm \n', smstd, smext)
                fprintf('Temporal smoothing kernel = %4.2f ms\n',Opts.tm*self.dt*1e3)
            end
        end
        
        function pad_axes(self)
            %%% pads the scan axes so the edge kernels will fit if the ROI
            %%% is on the edges of the image
            self.frame_size = self.frame_size + self.pad_inside + self.pad_outside;
            if strcmp(self.coord, 'polar')
                self.r = self.r(1)-self.dr*self.pad_inside(1): self.dr: ...
                    self.r(end)+self.dr*self.pad_outside(1);
                self.theta = self.theta(1)-self.dtheta*self.pad_inside(2): ...
                    self.dtheta: self.theta(end)+self.dtheta*self.pad_outside(2);
            elseif strcmp(self.coord, 'cartesian')
                self.z = self.z(1)-self.dz*self.pad_inside(1): self.dz: ...
                    self.z(end)+self.dz*self.pad_outside(1);
                self.x = self.x(1)-self.dx*self.pad_inside(2): self.dx: ...
                    self.x(end)+self.dx*self.pad_outside(2);
                self.origin = [self.z(1), self.x(1)];
            end
        end
        
        function crop_axes(self)
            %%% crop axes with ROI
            roi = self.ROI;  % to keep lines shorter
            if strcmp(self.coord, 'polar')
                self.r = self.r(roi(2):roi(2)+roi(4));
                self.theta = self.theta(roi(1):roi(1)+roi(3));
            elseif strcmp(self.coord, 'cartesian')
                self.z = self.z(roi(2):roi(2)+roi(4));
                self.x = self.x(roi(1):roi(1)+roi(3));
            end
            self.frame_size = [roi(3), roi(4)];
        end
        
        function calc_cartesian_axes(self, dx, dz)
            %%% creates cartesian axes for scan conversion
            %%% dx / dz - pixel resolution in mm
            self.dx = dx;
            self.dz = dz;
            r_lims = [self.r(1), self.r(end)];
            theta_lims = [self.theta(1), self.theta(end)];
            bbox_x = r_lims'*sin(theta_lims);
            bbox_z = r_lims;%'*cos(theta_lims);
            xlims = [min(bbox_x(:)), max(bbox_x(:))]; %min, max in mm
            zlims = [min(bbox_z(:)), max(bbox_z(:))]; %min, max in mm
            self.x = xlims(1):self.dx:xlims(2);  % x axis  [mm]
            self.z_0 = zlims(1):self.dz:zlims(2);  % z axis  (with reference to virtual focus) [mm]
            self.z = self.z_0 + self.origin(1);  % z axis [mm]
            self.z_0 = self.z_0(self.z > 0); % don't want to scan convert behind trans
            self.z = self.z(self.z > 0); % don't want to scan convert behind trans
            
            
        end
        
        function scan_convert_bmodes(self)
            %%%scan_convert - scan converts bmode in polar coords onto a cartesianb grid defined by
            %%%res and origin to grid.
            if isempty(self.bmodes_polar)
                warning('No polar bmode found, making one automatically from frames')
                self.create_bmode_by_coherent_compounding()
                self.ensemble_average_bmode()
            end
            [R, THETA] = ndgrid(self.r, self.theta);
            [Z, X] = ndgrid(self.z_0, self.x);
            self.bmodes_cart = scan_convert(self.bmodes_polar, R, THETA, Z, X);
        end
        
        function create_bmode_by_coherent_compounding(self)
            %%% coherently compounds image data and stores in a different
            %%% variable
            if strcmp(self.coord, 'polar')  % polar
                self.bmodes_polar = abs(squeeze(sum(self.frames,3)));
            else  % cartesian
                self.bmodes_cart = abs(squeeze(sum(self.frames,3)));
            end
        end
        
        function create_bmode_by_incoherent_compounding(self)
            %%% coherently compounds image data and stores in a different
            %%% variable
            if strcmp(self.coord, 'polar')  % polar
                self.bmodes_polar = squeeze(sum(abs(self.frames),3));
            else  % cartesian
                self.bmodes_cart = squeeze(sum(abs(self.frames),3));
            end
        end
        
        function ensemble_average_bmode(self)
            %%% ensemble averages bmode data for display behind vectors
            if strcmp(self.coord, 'polar')  % polar
                temp = self.bmodes_polar;
            else  % cartesian
                temp = self.bmodes_cart;
            end
            if Opts.n_ave > 1
                [nz, nx, ~] = size(temp);
                temp = reshape(temp,nz,nx,Opts.n_ave,[]);
                temp = squeeze(mean(temp,3));
                self.num_frames_bmode = size(temp,3);
            end
            if Opts.nm_ave > 1
                temp = movmean(temp,Opts.nm_ave,3);
            end
            if strcmp(self.coord, 'polar')  % polar
                self.bmodes_polar = temp;
            else  % cartesian
                self.bmodes_cart = temp;
            end
        end
        
        function [output] = get_bmode(self, bmode_type)
            %%% gets bmode depending on whether its polar or cartesian
            if nargin < 2
                if strcmp(self.coord, 'polar')
                    bmode_type = 'polar';
                else
                    bmode_type = 'cartesian';
                end
            end
            if strcmp(bmode_type,'polar')
                if isempty(self.bmodes_polar)
                    self.create_bmode_by_coherent_compounding()
                end
                output = self.bmodes_polar;
            elseif strcmp(bmode_type, 'cartesian')
                if isempty(self.bmodes_cart)
                    self.create_bmode_by_coherent_compounding()
                end
                output = self.bmodes_cart;
            else
                error('bmode type not recognised, use polar or cartesian')
            end
        end
        
        function [image1] = get_ref_image_stack(self, index)
            %%% gets the first stack of images (at time tn)
            image1 = (real(self.frames(:,:,:,Opts.n_ave*(index-1)+1:...
                index*Opts.n_ave+Opts.nm_ave-1)));
            if Opts.useGPU
                image1 = gpuArray(image1);
                
            end
            image1 = reshape(image1, size(image1,1), size(image1,2), []);
            
        end
        
        function [image2] = get_target_image_stack(self, index)
            %%% gets the second stack of images (at time tn+frame_skip)
            image2 = (real(self.frames(:,:,:,Opts.n_ave*(index-1)+Opts.frameSkip+1:...
                index*Opts.n_ave+Opts.frameSkip+Opts.nm_ave-1)));
            if Opts.useGPU
                image2 = gpuArray(image2);
            end
            image2 = reshape(image2, size(image2,1), size(image2,2), []);
        end
        
        %% mask methods
        function calcPadding(self)
            %%% calcs the necessary padding to the ROI and Frames/imMask so that
            %%% kernel centers fit on edges of ROI
            self.ROI = [self.ROI(1)-Opts.kSize(1,2)/2, self.ROI(2)-Opts.kSize(1,1)/2,...
                self.ROI(3)+Opts.kSize(1,2), self.ROI(4)+Opts.kSize(1,1)];
            self.pad_inside = flip(self.ROI(1:2));
            self.pad_inside(self.pad_inside > 0) = 0;
            self.pad_inside = abs(self.pad_inside);
            self.pad_outside = self.frame_size(1:2) - flip(self.ROI(3:4)+self.ROI(1:2));
            self.pad_outside(self.pad_outside > 0) = 0;
            self.pad_outside = abs(self.pad_outside);
            
            self.ROI = self.ROI + [flip(self.pad_inside)+1, -1, -1];
        end
        
        function pad_data(self)
            %%% padData Pads the mask and frames with values calculated from calc_padding.
            % pad frames
            self.frames = padarray(self.frames, double(self.pad_inside), 'symmetric', 'pre'); % pad
            self.frames = padarray(self.frames, double(self.pad_outside), 'symmetric', 'post'); % pad
            mask = self.get_mask();
            % pad mask
            mask = padarray(mask, double(self.pad_inside), 0, 'pre'); % pad
            mask = padarray(mask, double(self.pad_outside), 0, 'post'); % pad
            self.set_mask(mask);
            % now update axes
            self.pad_axes();
        end
        
        function crop_frames(self)
            %%% crops frames to ROI
            % crop frames
            self.frames = self.frames(self.ROI(2):self.ROI(2)+self.ROI(4),...
                self.ROI(1):self.ROI(1)+self.ROI(3),:,:);
            % crop mask
            mask = self.get_mask();
            mask = mask(self.ROI(2):self.ROI(2)+self.ROI(4),...
                self.ROI(1):self.ROI(1)+self.ROI(3));
            self.set_mask(mask);
            % crop axes
            self.crop_axes();
        end
        
        function create_mask(self)
            %%% convenience function for creating mask - choses appropriate method for
            %%% making a mask depending on data source
            if contains(Opts.datasource, 'SAVFI')
               self.create_SAVFI_mask(); 
            else
                self.create_manual_mask();
            end
        end
        
        function create_SAVFI_mask(self)
           %%% creates a mask for SAVFI data using the reference grid information
           flow_grid = get_SAVFI_grid(self.filepath);
           mask = zeros(self.frame_size);
           [Z,X] = ndgrid(self.z, self.x);
           in_z = (Z >= flow_grid.z(1)*1000) & (Z <= flow_grid.z(end)*1000);
           if contains(self.filepath,'straight_vessel_90deg')
               in_x = (X >= -5) & (X <= 5);
           elseif contains(self.filepath, 'straight_vessel_105deg')
               in_z = (Z >= flow_grid.z(1)*1000 + X.*sind(15)) & (Z <= flow_grid.z(end)*1000 + X.*sind(15));
               in_x = (X >= -5) & (X <= 5);            
           elseif contains(self.filepath, 'spinning_disk')
               in_z = (X.^2 + (Z - 25).^2) <= 8^2;  % inside circle (plus 0.5 mm)
               in_x = ones(size(X));  % hack to work with convention of other phantoms      
           else
                in_x = (X >= flow_grid.x(1)*1000) & (X <= flow_grid.x(end)*1000);
           end
           mask = in_z.*in_x;
           ROI = regionprops(mask, 'BoundingBox');
           self.ROI = fix(ROI.BoundingBox);
           mPoly = bwboundaries(mask);
           self.mPoly = mPoly{1};
           self.set_mask(mask);
        end
        
        function create_manual_mask(self)
            %%% manually draw a polygon defining the masked region
            if self.num_frames > 10
                filtframes = FIRFiltFrames(self.frames, [0.02, 1.2], self.cf, self.c, self.FR, 10);
                filtframes = filtframes(:,:,:,11:end);
            elseif self.num_frames > 9
                filtframes = ButterFiltFrames(self.frames, [0.02, 1.2], self.cf, self.c, self.FR, 3);
            else
                filtframes = self.frames - mean(self.frames,4);
            end
            filtframes = 20*log10(max(sum(abs(filtframes),3), [], 4));
            [mask, self.ROI, self.mPoly] = drawMask(filtframes);
            self.set_mask(mask);
            self.ROI = fix(self.ROI);
            self.ROI0 = self.ROI;  % keep copy of original ROI before modifying
            
        end
        
        function load_mask(self, dirname, filename)
            %%% load mask information from file
            maskpath = [dirname, '\', filename, '_mask.mat'];
            % load mask from file
            A = load(maskpath);
            self.set_mask(A.imMask);
            self.ROI = fix(A.ROI);
            self.ROI0 = self.ROI;  % keep copy of original ROI before modifying
            self.mPoly = A.mPoly;
        end
        
        function save_mask(self, dirname, filename)
            %%% save mask infomation to file
            maskpath = [dirname, '\', filename, '_mask.mat'];
            imMask = self.get_mask();
            ROI = self.ROI;
            mPoly = self.mPoly;
            save(maskpath,'imMask','ROI','mPoly');
        end
        
        function scan_convert_mask(self)
            %%%scan_convert_mask - scan converts mask in polar coords onto a cartesianb grid defined by
            %%%res and origin to grid.
            %%% need to input instance of ImageProperties
            assert(~isempty(self.imMask_polar), 'No polar mask found');
            [R, THETA] = ndgrid(self.r, self.theta);
            [Z, X] = ndgrid(self.z_0, self.x);
            
            temp = scan_convert(single(self.imMask_polar), R, THETA, Z, X);
            temp(isnan(temp)) = 0;
            self.imMask_cart = logical(temp);
        end
        
        function [mask] = get_mask(self)
            %%% gets mask whether polar or cartesian
            if strcmp(self.coord, 'polar')
                mask = self.imMask_polar;
            else
                mask = self.imMask_cart;
            end
        end
        
        function set_mask(self, mask)
            %%% gets mask whether polar or cartesian
            if strcmp(self.coord, 'polar')
                self.imMask_polar = mask;
            else
                self.imMask_cart = mask;
            end
        end
        
        function mask_frames(self)
            %%% convenience function to calculate and perform padding and
            %%% cropping of mask, frames and image axes.
            % first calculate if padding is required to fit largest kernel size
            % if ROI is close to edge of image data
            self.calcPadding();
            % pad frames and mask and axes
            self.pad_data();
            % crop frames and mask and axes to ROI
            self.crop_frames();
        end
        
        %% filtering methods
        function clutterFilter(self, mode)
            if nargin < 2
                mode = Opts.SlowTimeFilter;
            end
            switch mode
                case 'ButterIIR'
                    self.frames = ButterFiltFrames(self.frames, [Opts.min_vel, Opts.max_vel], self.cf, self.c, self.FR, fix(Opts.filterOrder/2));
                case 'FIR'
                    self.frames = FIRFiltFrames(self.frames, [Opts.min_vel, Opts.max_vel], self.cf, self.c, self.FR, Opts.filterOrder);      
                    % adjust frames to skip ramp-up time for filter
                    self.start_frame = self.start_frame+Opts.filterOrder; 
                    self.num_frames = self.num_frames - Opts.filterOrder;
                    self.frames = self.frames(:,:,:,Opts.filterOrder+1:end);
                    self.trim_frames() %update frame limits
                case 'SVDauto'
                    [self.frames, ~, ~, ~, ~] = svdfilterAuto(permute(self.frames, [1,2,4,3]), ...
                        Opts.tissue_thresh, Opts.rem_noise, Opts.noise_thresh);
                    self.frames = permute(self.frames, [1,2,4,3]); % put angles in dim 3 again
                case 'SVD'
                    [self.frames, ~, ~, ~] = svdfilter(permute(self.frames, [1,2,4,3]), ...
                        Opts.svd_min, Opts.svd_max);
                    self.frames = permute(self.frames, [1,2,4,3]); % put angles in dim 3 again
                case 'meanSub'
                    self.frames = self.frames - mean(self.frames,4);
                case 'movMeanSub'
                    self.frames = self.frames - movmean(self.frames,Opts.filterOrder,4);
                case 'None'
                otherwise
                    warning('Unknown filter specification, no slow-time filtering was performed')
            end
        end
        
        %% visulisation methods
        
        function play_bmode_with_mask(self, mode)
            %%% plays the bmode sequence with mask over it
            %%% mode - specifies if polar or cartesian sequence should be displayed
            %%% 'polar' or 'cartesian'
            if nargin < 2
                mode = 'Cartesian';
            end
            if strcmp(mode, 'polar')
                bmode = self.get_bmode('polar');
                mask = single(self.imMask_polar);
                hax = rad2deg(self.theta);
                vax = self.r;
                xtitle = 'Azimuth [deg]';
                ytitle = 'Depth [mm]';
            else
                bmode = self.get_bmode('cartesian');
                mask = single(self.imMask_cart);
                hax = self.x;
                vax = self.z;
                xtitle = 'X [mm]';
                ytitle = 'Z [mm]';
                
            end
            figure();
            hold on
            bmode = self.log_compression(bmode);
            fig = imagesc(hax,vax, bmode(:,:,1),[-40 0]);
            %plot mask
            mfig = imagesc(hax,vax, mask);
            alpha = 0.7.*mask;
            set(gca(), 'YDir', 'reverse')
            set(mfig, 'AlphaData', alpha)
            xlabel(xtitle);
            ylabel(ytitle);
            title(mode)
            %end plot mask
            colormap gray
            axis image
            for i = 1:size(bmode,3)
                fig.CData = bmode(:,:,i);
                pause(0.01)
            end
        end
        
        function play_frames(self, dt)
            if nargin < 2
                dt = 0.1; % pause time
            end
            %%% plays a preview of the frames
            if strcmp(self.coord, 'polar')
                bmode = self.get_bmode();
                hax = rad2deg(self.theta);
                vax = self.r;
                xtitle = 'Azimuth [deg]';
                ytitle = 'Depth [mm]';
                mode = 'Polar';
            else
                bmode = self.get_bmode();
                hax = self.x;
                vax = self.z;
                xtitle = 'X [mm]';
                ytitle = 'Z [mm]';
                mode = 'Cartesian';
            end
            figure();
            hold on
            bmode = self.log_compression(bmode);
            fig = imagesc(hax,vax, bmode(:,:,1),[-40 0]);
            colorbar()
            xlabel(xtitle);
            ylabel(ytitle);
            title(mode)
            set(gca,'YDir','reverse')
            %end plot mask
            colormap gray
            axis image
            for i = 1:size(bmode,3)
                fig.CData = bmode(:,:,i);
                pause(dt)
            end
        end
        
        function PIV_view_polar(self, pivfilt, drange, sf)
            %%% interactive visualisation of PIV vectors over bmode images
            %%% polar domain data
            %%% piv_filt - instance of PostProcess class
            %%% drange - dynamic range to show bmode
            if nargin < 4
                sf = 3.0;
            end
            if nargin < 3
                drange = [-50 0];
            end
            assert(~isempty(self.bmodes_polar), 'no polar bmode data');
            assert(pivfilt.is_spatial_disp, 'first convert polar piv data to spatial units');
            bmodes = self.log_compression(self.bmodes_polar(:,:,1:pivfilt.vdims(3)));
            X = rad2deg(pivfilt.Theta);
            Y = pivfilt.R;
            U = rad2deg(pivfilt.dTheta);
            V = pivfilt.dR;
            imlim_x = rad2deg([self.theta(1), self.theta(end)]);
            imlim_y = [self.r(1), self.r(end)];
            figure()
            vecscale(imlim_x, imlim_y,bmodes, drange, X, Y, U, V, sf)
        end
        
        function PIV_view_cart(self, pivfilt, drange, sf)
            %%% interactive visualisation of PIV vectors over bmode images
            %%% cartesian domain data
            %%% piv_filt - instance of PostProcess class
            %%% drange - dynamic range to show bmode
            if nargin < 4
                sf = 3.0;
            end
            if nargin < 3
                drange = [-50 0];
            end
            assert(~isempty(self.bmodes_cart), 'no cartesian bmode data');
            assert(pivfilt.is_spatial_disp || pivfilt.is_spatial_vel, 'first convert piv data to spatial units');
            bmodes = self.log_compression(self.bmodes_cart(:,:,1:pivfilt.vdims(3)));
            
            
            X = pivfilt.X;
            Y = pivfilt.Z;
            U = pivfilt.U;
            V = pivfilt.V;
            imlim_x = [self.x(1), self.x(end)];
            imlim_y = [self.z(1), self.z(end)];
            figure()
             vecscale(imlim_x, imlim_y,bmodes, drange, X, Y, U, V, sf)
%              vecscale(imlim_x, imlim_y,angles_dim, drange, X, Y, U, V, sf)
        end
        
        function save_bmode_avi(self, filepath, drange, coord)
            %%% saves bmode as avi
            %%% filepath - path to save avi to
            %%% drange - min and max intensity value (dB)
            %%% coord - 'polar' or 'cartesian'
            if nargin < 3
                drange = [-50 0];
                coord = self.coord;
            end
            bmode = self.get_bmode(coord);
            if strcmp(coord, 'polar')
                hax = rad2deg(self.theta);
                vax = self.r;
                xtitle = 'Azimuth [deg]';
                ytitle = 'Depth [mm]';
                mode = 'Polar';
            else
                hax = self.x;
                vax = self.z;
                xtitle = 'X [mm]';
                ytitle = 'Z [mm]';
                mode = 'Cartesian';
            end
            % Saving Bmodes
            v = VideoWriter(filepath, 'MPEG-4');
            v.FrameRate = 30;
            v.Quality = 100;
            open(v)
            fig = figure();
            hold on
            bmode = self.log_compression(bmode);
            im = imagesc(hax,vax, bmode(:,:,1),drange );
            set(gca,'YDir','reverse')
            xlabel(xtitle);
            ylabel(ytitle);
            title(mode)
            %end plot mask
            colormap gray
            axis image
            for i = 1:size(bmode,3)
                im.CData = bmode(:,:,i);
                frame = getframe(fig);
                writeVideo(v,frame); 
                title([mode '-' num2str(i)]);
            end
            close(v)
        end
        
        function [bmode_shape, bmode_type, bmode_path] = save_bmode_bin(self, filepath)
            %%% saves bmode images as binfil
            bmode = single(self.log_compression(self.bmodes_cart));
            bmode_shape = size(bmode);
            bmode_type = class(bmode);
            bmode_path = [filepath,'_bmode.bin'];
            fileID = fopen(bmode_path,'w');
            fwrite(fileID,bmode,bmode_type);
            fclose(fileID);
        end
        
    end %methods
end %class