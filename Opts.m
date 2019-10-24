classdef Opts
    properties (Constant)
        %% loading options
        datasource = 'Verasonics'; %'Verasonics','SAVFI', 'ustb'
        loadtype = 'matfiles';  % 'binfile', 'matfiles', 'uff'
        startFrame = 1;
        endFrame = 100;
        %% PIV options
        frameSkip = 1; %number of frames between frame comparisons for PIV (1=adjacent frames are compared)
%        ksize = single([[32,32];[32,32];[16,16];[16,16]]); % y,x kernel size without padding
       ksize = single([[16,16]]);%;[16,16]]); % increasing the number of kernels will not work
%         ksize = single([[32,32];[32,32];[16 16]]);
        padding = 1; % padding ksize to prevent wrap-around in fourier domain
        windowfunc = @boxcar;%  @blackman, @tukeywin, @boxcar, @hann % window to use for filtering kernels
        windowfunc_args = 1; %blackman/hann: 'periodic', 'symmetric'; %boxcar: 0-1; tukeywin: 0-1;  % additional arguments for window function if needed
        subPixFinder = 4; % 0 = none; 1 = gaussian 2x3 pixel; 2 = 2D (3x3) gaussian; 3 = centroid; 4 = 2x3 parabola; 5 = 2D (3x3) parabolic (dont use)
        n_ave = 4;  % number of frames used for correlation averaging
        nm_ave = 1; % number of frames used for correlation moving averaging
        imDeform = 'linear'; %image deformation method can be: linear, cubic, *spline - (only on CPU)
        stepDevisor = 4;  %devisor overlap for step size calculation. Must be multiple of 2 (or 1)
        step = Opts.ksize./Opts.stepDevisor;  % the step size in pixels.
        simOp = 'nxcc'; % 'nxcc', 'phasecorr' - type of similarity operator to use
        
        %% multi-pass options
        uni_args = {0.2, 4, 1};  %  universal outlier arguments - {epsilon, thresh, b} - see median_outlier.m
        mincorr = 0.05;%0.05; % minimum maximum normcorr value for a kernel
        useSmoothN = true;
        
         %% spatial weights
        SpatialWeights = 1; % vessel orientation restriction
        w_type = 'Huber';% 'Cauchy' 'Huber' 'Bisquare' 
        
        %% directional weights
        dc = 0; % apply directional weights
        dsize = 5;%15; % window size to smoothen the first order derivatives
        
        %% smoothing options
        smoothing = 2;  % do spatial smoothing during post processing
        sm = [0.5,0.5,3,3].*Opts.smoothing;  % smoothing kernel sizes
        tm = 3;  % temporal smoothing length
        remove_outliers = true;%true;
        use_universal_outlier = true;%true;
        fill_outliers = true;%true;
        useSmoothN_robust = true;%true;
        
        %% script options
        useGPU = 0;  % option to use GPU or not
        piv_save = 1;  % if we are going to save PIV data
        PIVVIEW = 0;  % if we want to view PIV results in Matlab
        live_view = 1; %show quiver plots as they are calculated (slow down calc)
        corr_over_quiv = 1;  % show quiver plot over correlation map each frame
        maskvecs = 2;  % number of vectors to mask for live view visualisation - 1 = no masking, 2 = mask every second
        piv_makeMovie = 1;  % make a movie using python script
        bmodeType = 'filt'; % 'filt' (show filtered data for bmode) 'unfilt' (unfiltered)
        debug = 0;
        
        %% tissue suppression options
        SlowTimeFilter = 'ButterIIR'; % 'ButterIIR' 'FIR' 'None' 'SVDauto', 'movMeanSub','meanSub', 'SVD'
        filterOrder = 12;  % filter order for ButterIIR, FIR and moving mean ensemble size for movMeanSub
        % frequency based filter specs
        min_vel = 0.04 % minimum velocity in m/s
        max_vel = 2; % maximum velocity in m/s
        % SVDauto filter specs
        tissue_thresh = 0.99;  % slope of singular value curve to cut-off low-order modes
        rem_noise = 0;  % whether we want to remove high order singular values
        noise_thresh = 0.9; % noise threshold if applicable
        % SVD filter cutoffs
        svd_min = 1;
        svd_max = 250;
        
        
        %% calculated options
        nIterations = sum(Opts.ksize(:,1)~=0); %window size and iterations 
        fin_size = Opts.ksize(end,:);
        procFrames = Opts.endFrame-Opts.startFrame+1; %number of frames to be processed. if 0 process all frames
        kSize = Opts.ksize.*Opts.padding;  %kernel size after padding
        %logcomp = 0;
        
    end
end
        