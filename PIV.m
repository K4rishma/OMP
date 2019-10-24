classdef PIV < handle
    % PIV class containing all PIV data, methods and properties
    properties (Access = public)
        image1          % reference image stack (vert, horz, time)
        image2          % target image stack (vert, horz, time)
        piv_start       % start frame for PIV
        piv_stop        % end frame for PIV
        num_vframes     % number of PIV frames
        minis           % indices (x, y, t) of first kernel
        maxis           % indices (x, y, t) of last kernel
        n_els           % number of kernels that fit in each direction (x, y, t)
        ss1s            % linear index matrix for vectorized operation (frame 1)
        ss2s            % linear index matrix for vectorized operation (frame 2)
        fdims           % size of a frame
        X               % center positions of kernels along horizontal axis (pixels)
        Z               % center positions of kernels along vertical axis (pixels)
        U               % displacement values in horizontal direction (pixels/frame)
        V               % displacement val ues in vertical direction (pixels/frame)
        vflags          % vector flag status (int8)
        corrMap         % maximum correlation per kernel
        findex          % current frame index (must be a better way?)
        Wf              % window filter for minimizing spectral overlap
        image_wfunc     % kernel window function object (1d)
        corrfunc        % correlation function - generic handle for assignment during piv_init
        multipass       % current index in multipass
        passes          % number of passes to make
        sw              % spatial weight
        %visualisation
        show_corr_live  % option for showing quiver over correlation map on each frame iteration
        lvcorr_fig      % figure handle for quiver over correlation map figure
        lvcorr_im       % image handle for correlation data for quiver over correlation map figure
        lvcorr_quiv     % handle to quiver data for quiver over correlation map figure
        lvcorr_sf       % scaling factor for quivers
        lvcorr_mf       % mask factor so we can show less vectors (reduce density)
        log_compression % quick function for taking envelope and log compressing
        show_quiv_live  % option for showing quiver over bmode image
        lvquiv_fig      % figure handle for quiver over image figure
        lvquiv_im       % image handle for image which quiver is overlaid
        lvquiv_sf       % scaling factor for quivers
        lvquiv_quiv     % handle to quiver data for quiver over bmode figure
        lvquiv_mf       % mask factor so we can show less vectors (reduce density)    
    end
    %% constructor
    methods (Access = public)
        function self=PIV()
            self.minis = {};
            self.maxis = {};
            self.n_els = {};
            self.ss1s = {};
            self.ss2s = {};
            self.Wf = {};
            self.log_compression = @(x) 20*log10(abs(x/max(x(:))));  % handy function for log compression    
        end
    end
    %% methods
    methods
        function init(self, tend, fdims)
            %Initializes the variables for PIV that only need to be computed once and
            %not on every loop
            self.piv_start = 1;
            self.piv_stop = floor((tend-Opts.frameSkip-(Opts.nm_ave-1))/(Opts.n_ave))-1;
            self.num_vframes = self.piv_stop-self.piv_start+1;
            self.fdims = fdims;
            self.make_window_function()
            % assign correlation function
            switch Opts.simOp
                case 'nxcc'
                    self.corrfunc = @normxcorr;
                case 'phasecorr'
                    self.corrfunc = @phasecorr;
            end
            
            % vectorize block matching process for each window size
            % iteration
            for i = 1:size(Opts.kSize,1)
                cas = [];%[];
                [mini, maxi, n_el] = get_extents([Opts.kSize(i,:), 1], [Opts.step(i,:), 1], self.fdims,cas);
                self.minis{i} = mini; self.maxis{i} = maxi; self.n_els{i} = n_el;
                self.ss1s{i} = kernalize(mini, maxi, [Opts.step(i,:), 1], n_el, self.fdims, [Opts.kSize(i,:), self.fdims(3)]);
                
                cas = 2;%2
                [mini2, maxi2, n_el2] = get_extents([Opts.kSize(i,:).*2, 1], [Opts.step(i,:), 1], self.fdims, cas);
                %self.minis{i} = mini; self.maxis{i} = maxi; self.n_els{i} = n_el2;
                self.ss2s{i} = kernalize(mini2, maxi2, [Opts.step(i,:), 1], n_el2, self.fdims, [Opts.kSize(i,:).*2, self.fdims(3)]);
                
                % make window function for this iteration
                Wfz = self.image_wfunc(Opts.kSize(i,1));
                Wfx = self.image_wfunc(Opts.kSize(i,2));
                self.Wf{i} = Wfz*Wfx';
                if i > 1  % target window has slightly different vectorization procedure
                    %- only needed for 2nd iteration onwards
                    dims = maxi(1:2) + Opts.kSize(i,:) - mini(1:2);
                    self.ss2s{i-1} = kernalize([1,1,1], [Opts.step(i,:).*n_el(1:2), self.fdims(3)], ...
                        [Opts.step(i,:), 1], n_el, [dims, self.fdims(3)], [Opts.kSize(i,:), self.fdims(3)]);
                end
            end
            % initialize arrays
            x = [mini(2):Opts.step(end,2):maxi(2)]+Opts.kSize(end,2)/2;
            z = [mini(1):Opts.step(end,1):maxi(1)]+Opts.kSize(end,1)/2;
            [self.X, self.Z] = meshgrid(x,z);
            self.U = zeros(size(self.X,1), size(self.X,2), self.num_vframes, 'single');
            self.V = zeros(size(self.U), 'single');
            self.vflags = ones(size(self.U), 'int8');
            self.corrMap = zeros(size(self.U), 'single');
            % initialize visualisations if needed
            if Opts.corr_over_quiv
                sf = 2;
                self.live_view_corr_init(sf,Opts.maskvecs)
            end
            if Opts.live_view
                sf = 3;
                self.live_view_quiv_init(sf, [-60 0],Opts.maskvecs)
            end
            
        end
        
        function make_window_function(self)
            %%% creates the window function defined by user along one dimension.
            %%% Needs to be called once per dimension and multiplied to
            %%% produce 2D window function
            %%% ksize = size of kernel in pixels along one spatial dimenion
            
            % locally assign variable names for readibility / writeability
            wf_args = Opts.windowfunc_args;
            wf_str = func2str(Opts.windowfunc);
            function win1d = boxcar_custom(ksize)
                assert(isnumeric(wf_args), 'Opts.windowfunc_args should be numeric for boxcar/rectwin window function');
                assert((wf_args <= 1) && (wf_args >=0), 'Opts.windowfunc_args should be between 0 and 1 for boxcar/rectwin window function');
                win1d = zeros(ksize,1);
                idx0 = ceil((ksize/2)-(ksize*wf_args/2))+1; % start indices for top-hat
                idx1 = floor((ksize/2)+(ksize*wf_args/2)); % end indices for top-hat
                win1d(idx0:idx1) = 1; % set top hat
            end
            function win1d = blackman_custom(ksize)
                assert((strcmp(wf_args, 'periodic') || strcmp(wf_args, 'symmetric')), ...
                    'Opts.windowfunc_args should be a string ''periodic'' or ''symmetric'' for blackman window type');
                win1d = blackman(ksize, wf_args);
            end
            function win1d = tukeywin_custom(ksize)
                assert(isnumeric(wf_args), 'Opts.windowfunc_args should be numeric for tukeywin window function');
                assert((wf_args <= 1) && (wf_args >=0), 'Opts.windowfunc_args should be between 0 and 1 for tukeywin window function');
                win1d = tukeywin(ksize, wf_args);
            end
            function win1d = hann_custom(ksize)
                assert((strcmp(wf_args, 'periodic') || strcmp(wf_args, 'symmetric')), ...
                    'Opts.windowfunc_args should be a string ''periodic'' or ''symmetric'' for hann window type');
                win1d = hann(ksize, wf_args);
            end
            switch wf_str  % switch between known matlab window functions (not all implemented yet)
                % this is needed because I don't know how to pass a
                % variable number of arguments to a generic function handle
                % - any suggestions appreciated
                case 'rectwin'
                    if ~isempty(wf_args)
                        self.image_wfunc = @boxcar_custom;
                    else
                        self.image_wfunc = @rectwin;
                    end
                case 'boxcar'
                    if ~isempty(wf_args)
                        self.image_wfunc = @boxcar_custom;
                    else
                        self.image_wfunc = @rectwin;
                    end
                case 'blackman'
                    if ~isempty(wf_args)
                        self.image_wfunc = @blackman_custom;
                    else
                        self.image_wfunc = @blackman;
                    end
                case 'tukeywin'
                    if ~isempty(wf_args)
                        self.image_wfunc = @tukeywin_custom;
                    else
                        self.image_wfunc = @tukeywin;
                    end
                case 'hann'
                    if ~isempty(wf_args)
                        self.image_wfunc = @hann_custom;
                    else
                        self.image_wfunc = @hann;
                    end
                otherwise
                    self.image_wfunc = Opts.windowfunc;
            end
        end
        
        function process(self,scan, findex)
            % init stuff
            self.findex = findex;
            self.passes = size(Opts.kSize,1);
            for passnum= 1:self.passes
                self.multipass = passnum;
                % init variables per pass
                step = Opts.step(self.multipass,:);
                ksize = Opts.kSize(self.multipass,:);
                
                mask = scan.get_mask();
                mask(mask>1)=1;
                
                mini = self.minis{self.multipass};
                maxi = self.maxis{self.multipass};
                n_el = self.n_els{self.multipass};
                ss1 = self.ss1s{self.multipass};
                winFilt = self.Wf{self.multipass};
                typevector= ones(n_el(1),n_el(2), 'int8');
                
                if self.multipass > 1 % regularize previous passes (skipped for first pass)
                    %init vectorization for target kernel
                    ss2 = self.ss2s{self.multipass-1};  % one less ss2 (only starts on second iteration)
                    
                    % now regularize before deforming image
                    if Opts.debug,figure(1); imagesc(self.calc_mag(utable,vtable)); colorbar;end
                    [utable, vtable] = self.regularize_multipass(utable, vtable, Wcorr, typevector);  
                    if Opts.debug,figure(2); imagesc(self.calc_mag(utable,vtable)); colorbar;end
                    
                    % deform target image
                    [image2_deformed, utable, vtable] = ImageDeform(self.image2, xtable, ytable, utable, ...
                        vtable,n_el, mini, maxi, step, ksize);  % deform image2
                else
                    %init vectorization for target kernel - on first
                    %iteration same as reference kernel
                    ss2 = self.ss2s{self.multipass};
                    % init velocities to be 0
                    utable = zeros(n_el(1:2));
                    vtable = zeros(n_el(1:2));
                    % reference target image - first iteration no image
                    % deformation
                    image2_deformed = self.image2;
                end
                
                % divide images by small pictures
                image1_cut = self.image1(ss1).*winFilt; % create reference kernel stack
%                 image2_cut = image2_deformed(ss2);%.*winFilt; % create target kernel stack
                image2_cut = self.image2(ss2);%.*winFilt;
                
                % get correlation maps
                corrmaps = squeeze(self.corrfunc(image1_cut,image2_cut));
                
                % applying directional weights
                if Opts.dc == 1
                    acquisitions = scan.frames;
                    angles = anglesFromDerivatives(acquisitions);
                    corrmaps = directional_constraint(ss2,corrmaps,angles);
                    % correct the multiplication
                end
                %apply mask
                ii = ~mask(squeeze(ss1(round(ksize(1)/2+1), round(ksize(2)/2+1), 1, :))); % image pixel indices of mask
                jj = ~mask((mini(1):step(1):maxi(1))+round(ksize(1)/2), (mini(2):step(2):maxi(2))+round(ksize(2)/2)); % vector pixel indices of mask
                corrmaps(:,:, ii) = NaN; % mask correlation values
                typevector(jj) = 0; % mask vectors
                
                % subpixel fit
                %[vector1, maxVals,hort1,vert1] = subpix_fit(corrmaps, ksize);
                %vector1(ii,:)  = NaN;
                
%                  kk = repmat(~ii,1,2);
%                 hort1 = (hort1(~ii)).'; vert1 = (vert1(~ii)).';
%                 B = vector1(kk);
%                 B = reshape(B, [], 2);
%                vector1(ii,:)  = NaN; vector = vector1;
               
                 meanIm1 = mean(mean(mean(image1_cut)));
                 meanIm2 = mean(mean(mean(image2_cut)));
                meanIm1 = 0;
                meanIm2 = 0;
                image1m = squeeze(mean(image1_cut, 3));
    
                
                image2m = squeeze(mean(image2_cut , 3));
                image1m = image1m - meanIm1;
                image2m = image2m - meanIm2;
                
                if findex == 1
                
                Ya = repmat([ksize(1)/2: -1: -(ksize(1)/2 )]', 1, ksize(2)+1);
                Xa = repmat([-ksize(2)/2: 1: (ksize(2)/2)], ksize(1)+1, 1);
                
                angular_pos = atan2d(Ya,Xa);
                
                acquisitions = scan.frames;
                angles = anglesFromDerivatives(acquisitions);
%                 figure; imagesc(angles);
                
                w_type = 'Huber';
%                 s = 1/4;
                switch w_type
                       
                    case 'Cauchy'
                        s = 1/16;
                        cauchy = @(x) min(2/(1 + (x*s)^2),1);
                        [sw,med_angles] = spatialDistribution(angles, cauchy, angular_pos,ss1);
                    case 'Huber'
                        p = 1.5; s = 1/16;% p = 1.5
                        huber = @(x) min(p/abs((x*s)),1);
                        [sw,med_angles] = spatialDistribution(angles, huber, angular_pos,ss1);
                    case 'Bisquare'
                        c = 22; s = 1/16;
                        bisquare = @(x) max((1 - ((x*s)/c)^2)^2,0);
                        [sw,med_angles] = spatialDistribution(angles, bisquare, angular_pos,ss1);
                    otherwise
                        disp('unknown error');
                end
                
                self.sw = sw;
                
                else
                    
                    sw = self.sw;
                end
% angles1 = reshape(med_angles.', n_el(1), n_el(2));
% figure; imagesc(angles1);
% hax = [scan.x(1), scan.x(end)];
% vax = [scan.z(1), scan.z(end)];
% figure();
% hold on;
% fig = imagesc(hax,vax, real(scan.frames(:,:,1,1)));
% %plot mask
% %mfig = imagesc(hax,vax, mask);
% alpha = 0.7.*mask;
% set(gca(), 'YDir', 'reverse')
%  %set(mfig, 'AlphaData', alpha)
% 
% %end plot mask
% dim1 = size(scan.frames,1);
% dim2 = size(scan.frames,2);
% angle_dim = NaN(dim1,dim2);
% 
% for indsz = 1:1:size(self.Z,1)
% for indsx = 1:1:size(self.X,2)
%  
%     
% first = self.Z(indsz,indsx) : self.Z(indsz,indsx) + Opts.step(1) - 1;
% second = self.X(indsz,indsx) : self.X(indsz,indsx) + Opts.step(2) - 1;
% angle_dim(first, second) = angles1(indsz,indsx);
% end
% end
% 
% mask_v = double(mask);
% mask_v(mask==0) = NaN;
% fig.CData = angle_dim.*mask_v;%angle_dim;
% 
%                 hold off;
%             sw(:,:,:) = 1;%temporary
%                 [vector,hort,vert] = RobustOMP(image1m, image2m, sw, ii);
%                 kk = repmat(~ii,1,2);
%                 hort = hort(~ii); vert = vert(~ii);
%                 A = vector(kk);
%                 A = (reshape(A, [],2));
%                 [vector,hort,vert] = correlation1(image1_cut, image2_cut, sw, ii); 
                [vector, max_value] = correlation2(image1_cut, self.image2(ss2), sw, ii); % ss2 is different from ss1
                % new base tables
                xtable = repmat((mini(2):step(2):maxi(2))+ksize(2)/2, length(mini(1):step(1):maxi(1)), 1);
                ytable = repmat(((mini(1):step(1):maxi(1))+ksize(1)/2)', 1, length(mini(2):step(2):maxi(2)));
                
                if Opts.debug
                    figure(10);
                    sf = 5;
                    if (rem(ksize,2) == 0) %for the subpixel displacement measurement
                        SubPixOffset=1;
                    else
                        SubPixOffset=0.5;
                    end
                    debugPIVMulti(self.image1,image1_cut,image2_cut,corrmaps,vector, ...
                        xtable, ytable, SubPixOffset, [-50 0],sf)
                end
                
                vector = reshape(vector, [size(xtable) 2]);
                %Wcorr = reshape(maxVals, size(xtable));
                Wcorr = reshape(max_value, size(xtable));
                Wcorr(isnan(Wcorr)) = 0;
                Wcorr(Wcorr<0) = 0;
                 utable(Wcorr<Opts.mincorr) = nan;
                 vtable(Wcorr<Opts.mincorr) = nan;
                
                utable = utable+vector(:,:,1);
                vtable = vtable+vector(:,:,2);
                
            end  %multipass
            
            self.U(:,:,self.findex) = gather(utable);
            self.V(:,:,self.findex) = gather(vtable);
            self.vflags(:,:,self.findex) = gather(typevector);
            self.corrMap(:,:,self.findex) = gather(Wcorr);
            
        end
        
        function [utable, vtable] = regularize_multipass(self, utable, vtable, Wcorr, typevector)
            %%% regularisation of previous iterations displacements before
            %%% deforming target image
            
            % median test - look here if high velocities go missing
            args = Opts.uni_args;  % for some reason matlab doesnt let you use Opts.uni_args directly in function call?
            [utable, vtable] = median_outlier(utable, vtable, args{:});
            %replace nans
            utable(typevector==0) = 0;
            vtable(typevector==0) = 0;
%             utable(Wcorr<Opts.mincorr) = nan;
%             vtable(Wcorr<Opts.mincorr) = nan;
            %                 Wcorr(typevector<=0) = 1;
            if Opts.useGPU
                utable=gpuArray(inpaint_nans(gather(utable),4));
                vtable=gpuArray(inpaint_nans(gather(vtable),4));
            else
                utable=single(inpaint_nans(double(utable),4));
                vtable=single(inpaint_nans(double(vtable),4));
            end
            %smooth predictor - look here for improvements
            if Opts.useSmoothN
                try
                    if self.multipass<self.passes-1
                        utable = smoothn(utable, 0.1); %stronger smoothing for first passes
                        vtable = smoothn(vtable, 0.1);
                    else
                        utable = smoothn(utable,0.05); %weaker smoothing for last pass
                        vtable = smoothn(vtable,0.05);
                    end
                catch
                    error('could not find smoothn, see: Garcia D, http://www.biomecardio.com/pageshtm/publi/csda10.pdf')
                end
            else
                % gaussian kernel
                filt_handle=fspecial('gaussian',3,1);
                utable=imfilter(utable,filt_handle,'replicate');
                vtable=imfilter(vtable,filt_handle,'replicate');
            end
        end
        
        
        %% visualisation methods
        function live_view_corr_init(self, scale_factor, mask_factor)
            %%% live_view_quiv Plots quiver for a frame over the mean of the image stack
            %%% used to calculate the quiver
            if nargin < 3
                mask_factor = 1;
            end
            self.lvcorr_mf = mask_factor;
            self.lvcorr_sf = scale_factor;
            self.lvcorr_fig = figure();
            self.lvcorr_im = imagesc(self.corrMap(:,:,1), [0,1]);
            colormap hot
            colorbar
            hold on
            indsz = 1:self.lvcorr_mf:size(self.X,1);
            indsx = 1:self.lvcorr_mf:size(self.X,2);
            self.lvcorr_quiv = quiver(indsx, indsz,...
                self.U(indsz, indsx,1).*self.lvcorr_sf, self.V(indsz, indsx,1).*self.lvcorr_sf,0);
            title(1)
            drawnow;
            
        end
        
        function live_view_corr_update(self)
            %%% updates data in quiver over correlation plot
            self.lvcorr_im.CData = self.corrMap(:,:,self.findex);
            indsz = 1:self.lvcorr_mf:size(self.X,1);
            indsx = 1:self.lvcorr_mf:size(self.X,2);
            self.lvcorr_quiv.UData = self.U(indsz, indsx,self.findex).*self.lvcorr_sf;
            self.lvcorr_quiv.VData = self.V(indsz, indsx,self.findex).*self.lvcorr_sf;
            title(self.findex)
            drawnow;
        end
        
        function live_view_quiv_init(self, scale_factor, drange, mask_factor)
            %%% Plots quiver for a frame over the mean of the image stack
            %%% %used to calculate the quiver
            if nargin < 4
                mask_factor = 1;
            end
            self.lvquiv_mf = mask_factor;
            self.lvquiv_sf = scale_factor;
            self.lvquiv_fig = figure();
            im = ones(self.fdims(1:2));
            self.lvquiv_im = imagesc(self.log_compression(im), drange);
            colormap gray
            hold on
            indsz = 1:self.lvquiv_mf:size(self.X,1);
            indsx = 1:self.lvquiv_mf:size(self.X,2);
            self.lvquiv_quiv = quiver(self.X(indsz,indsx),self.Z(indsz,indsx),...
                self.U(indsz,indsx,1).*self.lvquiv_sf, self.V(indsz,indsx,1).*self.lvquiv_sf,0);
            title(1)
            drawnow;
        end
        
        function live_view_quiv_update(self)
            %%% updates data in quiver over image figure - needs an image
            %%% input
            if Opts.useGPU
                imdata = gather(self.log_compression(mean(self.image1,3)));
            else
                imdata = self.log_compression(mean(self.image1,3));
            end
            self.lvquiv_im.CData = imdata;
            indsz = 1:self.lvquiv_mf:size(self.X,1);
            indsx = 1:self.lvquiv_mf:size(self.X,2);
            self.lvquiv_quiv.UData = self.U(indsz,indsx,self.findex).*self.lvquiv_sf;
            self.lvquiv_quiv.VData = self.V(indsz,indsx,self.findex).*self.lvquiv_sf;
            title(self.findex)
            drawnow;
        end
        
        function clean_up(self)
            %%% closes live view windows if they are open
            if ~isempty(self.lvcorr_fig)
                close(self.lvcorr_fig)
            end
            if ~isempty(self.lvquiv_fig)
                close(self.lvquiv_fig)
            end
        end
        
        function [mag] = calc_mag(~, u, v)
           %%% calculates magnitude of displacement vectors
           mag = sqrt(u.^2 + v.^2);
        end
        
    end
    
    
end

