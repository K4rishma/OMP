clearvars
clc
doBatch = false;
doBatchallPhrase = false;

if doBatchallPhrase
    rootdir = uigetdir('K:\Results\Phrase\', 'Pick top level folder, all 1st level subfolders will be processed');
    subdirs = dir([rootdir '\Patient*']);  % this might be different for you (my folders started with 2018
    subdirs = {subdirs([subdirs.isdir]).name}; % make sure folders are folders
    for sdidx = 1:size(subdirs,2)  % loop through all subfolders (1 level lower)
        subsubdirs = dir([rootdir '\' subdirs{sdidx} '\2018*']);
        subsubdirs = {subsubdirs([subsubdirs.isdir]).name}; % make sure folders are folders
        for ssdidx = 1:size(subsubdirs,2)  % loop through all subfolders (1 level lower)
            DATA_PATH = [rootdir '\' subdirs{sdidx} '\' subsubdirs{ssdidx} '\'];
            temp = dir([DATA_PATH '*_1_info.mat']);
            FILENAME = temp(1).name;
            RunPIVFunction(FILENAME, DATA_PATH)
        end
    end
    
elseif doBatch
    % get root directory for batch processing
    rootdir = uigetdir('K:\Results\Phrase\', 'Pick top level folder, all 1st level subfolders will be processed');
    subdirs = dir([rootdir '\2018*']);  % this might be different for you (my folders started with 2018
    subdirs = {subdirs([subdirs.isdir]).name}; % make sure folders are folders
    for sdidx = 1:size(subdirs,2)  % loop through all subfolders (1 level lower)
        if strcmp(Opts.loadtype, 'matfiles')  % if we are loading matfiles
            DATA_PATH = [rootdir '\' subdirs{sdidx} '\'];
            temp = dir([DATA_PATH '*_1_info.mat']);
        elseif strcmp(Opts.loadtype, 'binfile') % if we are loading binfiles
            DATA_PATH = [rootdir '\' subdirs{sdidx} '\conc\'];
            temp = dir([DATA_PATH '*.bin']);
        end
        FILENAME = temp(1).name;
        RunPIVFunction(FILENAME, DATA_PATH)
    end
else
    if strcmp(Opts.datasource, 'Verasonics')
        if strcmp(Opts.loadtype, 'matfiles')
            [FILENAME,DATA_PATH,~] = uigetfile({'*_1_info.mat',  'All';}, ...
                  'Pick a file','K:\Results\Phrase\');
        elseif strcmp(Opts.loadtype, 'binfile')
            [FILENAME,DATA_PATH,~] = uigetfile({'*.bin',  'All';}, ...
                  'Pick a file','K:\Results\Phrase\');
        end 
    elseif strcmp(Opts.datasource, 'SAVFI')
        if strcmp(Opts.loadtype, 'matfiles')
            [FILENAME,DATA_PATH,~] = uigetfile({'*.mat',  'MATLAB file';}, ...
                      'Pick a file','K:\Results\SAVFI\SAVFI-data\');
        elseif strcmp(Opts.loadtype,'uff')
            [FILENAME,DATA_PATH,~] = uigetfile({'*.uff',  'Ultrasound File Format file';}, ...
                      'Pick a file','K:\Results\SAVFI\SAVFI-data\');
        end
    end
    RunPIVFunction(FILENAME, DATA_PATH)
    
end
%%
function RunPIVFunction(FILENAME, DATA_PATH)
%%  Load Data
    clc
    addpath('statsLib');  % this assumes 2DEchoPIV folder is working dir.
    addpath('Plotting');  % this assumes 2DEchoPIV folder is working dir.
    %load data - can change to your own. Frames must be MxNxTime
    filepath = fullfile (DATA_PATH, FILENAME);
    scan = ImageProperties();
    if strcmp(Opts.datasource,'Verasonics')  % load verasonics files
        if strcmp(Opts.loadtype, 'binfile')  
            scan.loadFromBinFileWrapper(filepath);
        elseif strcmp(Opts.loadtype, 'matfiles')
            scan.loadFromMatfilesWrapper(filepath);
        end    
    elseif strcmp(Opts.datasource, 'SAVFI')  % load SAVFI data
        if strcmp(Opts.loadtype, 'matfiles')
            scan.loadSAVFIdata(filepath)
        elseif strcmp(Opts.loadtype, 'uff')  % load SAVFI data from uff file
            scan.loadSAVFIustbdata(filepath)
        end
    else
        error('Unknown data source, please program your own');
    end

%     scan.choose_angles(1:3);
%     scan.frames = sum(scan.frames,3);
%     scan.num_ang = 1;
    
    % print some info relating PIV options to data
    scan.printPIVInfo   

    % Select/Load and Apply Mask 
    try
        assert(~contains(Opts.datasource, 'SAVFI'))
        scan.load_mask(DATA_PATH, FILENAME(1:find(FILENAME=='_')-1));
        % apply the mask
        scan.mask_frames();
    catch
        scan.create_mask();
        scan.save_mask(DATA_PATH, FILENAME(1:find(FILENAME=='_')-1));
        % apply the mask
        scan.mask_frames();
    end

    % if we want unfiltered bmode then we should create it before filtering
    if strcmp(Opts.bmodeType, 'unfilt')  
        % coherent compounding for bmode background
        scan.create_bmode_by_coherent_compounding();  
        % ensemble average bmode background by same amount as vector corr compounding
        scan.ensemble_average_bmode(); 
    end
  
    %  Filter Data
    scan.clutterFilter();
    
    if strcmp(Opts.bmodeType, 'filt')
        scan.create_bmode_by_coherent_compounding();  %coherent compounding for bmode background
        scan.ensemble_average_bmode(); % ensemble average bmode background by same amount as vector corr compounding
    end
    
    % Init PIV
    piv = PIV();
    
    fdims = [size(scan.frames(:,:,1,1)), Opts.n_ave*scan.num_ang+Opts.nm_ave-1];
    piv.init(scan.num_frames, fdims);

    % Process PIV
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
   
    % Scan convert bmode images
    if strcmp(scan.coord,'polar')
        scan.calc_cartesian_axes(scan.dr, scan.dr);
        scan.scan_convert_bmodes(); % convert frames
        scan.scan_convert_mask(); % convert mask   
    end

    % Post processing
    pivfilt = PostProcess(piv, scan);
    pivfilt.postprocess();
    % debug
%     scan.play_bmode_with_mask('cartesian')
    % scan covert vector maps
    if strcmp(scan.coord, 'polar')
        % convert vectors to cartesian (still on polar grid)
        pivfilt.scan_convert_vectors();
        % interpolate onto cartesian grid
        pivfilt.scan_convert_grid();
    end
    pivfilt.convert_disp_to_vel();
    % debug
%     pivfilt.play_velmag();
    % View PIV Results - Polar
    if Opts.PIVVIEW
       scan.PIV_view_polar(pivfilt, [-60 0], 4.0);
       scan.PIV_view_cart(pivfilt, [-60 0], 4.0);
    end

    % Save data for python plotting
    if Opts.piv_save
        % determine save path
        descstr = '_k';
        for k=1:Opts.nIterations
            descstr = [descstr ,num2str(Opts.ksize(k,1)),'x'];
        end
        descstr = descstr(1:end-1);
        descstr = [descstr '_pd' num2str(Opts.padding) func2str(Opts.windowfunc) num2str(Opts.windowfunc_args)];
        if contains(Opts.SlowTimeFilter, 'ButterIIR') || contains(Opts.SlowTimeFilter, 'FIR')
            descstr = [descstr Opts.SlowTimeFilter 'O' num2str(Opts.filterOrder) 'mv' num2str(Opts.min_vel)];
        elseif contains(Opts.SlowTimeFilter, 'SVD')
            descstr = [descstr Opts.SlowTimeFilter 'min' num2str(Opts.svd_min) 'max' num2str(scan.num_frames-Opts.svd_max)];
        else
            descstr = [descstr Opts.SlowTimeFilter];
        end
        descstr = [descstr '_cave' num2str(Opts.n_ave) '_cmave' num2str(Opts.nm_ave) '_vave', num2str(Opts.tm)];  % cave=correlation averaging, cnave = correlation moving averaging, vave = vector moving averaging
        descstr = [descstr Opts.simOp];
        descstr(strfind(descstr,'.')) = ',';
        pathparts = split(DATA_PATH, '\');
        pathpartskeep = {pathparts{1:3}, 'PIV', pathparts{4:length(pathparts)-1}};
        RESULTS_PATH = join(pathpartskeep,'\');
        RESULTS_PATH = RESULTS_PATH{1};
        mkdir(RESULTS_PATH)
        [~, fname, ~] =fileparts(fullfile(DATA_PATH, FILENAME));
        savepath = [RESULTS_PATH,'\', fname, descstr];
        % save
        pivfilt.save_piv_data(savepath, scan);
        % save SAFVI data if applicable
        if strcmp(Opts.datasource, 'SAVFI')
            n_emissions = Opts.nm_ave*Opts.n_ave*scan.num_ang*Opts.tm;
            if contains(savepath, 'carotid') && contains(savepath,'testing')
                pivfilt.save_SAVFI_results(savepath, n_emissions)
            else
                pivfilt.get_SAVFI_results(savepath, n_emissions)
            end
        end
        fprintf('Finished Saving')
    end
    %
    if Opts.piv_makeMovie
        systemCommand = ['python Python\PathlinePlot.py ', '"', savepath, '"',' ', fname];
        system(systemCommand)
%         systemCommand = ['python Python\QuiverPlot.py ', '"', savepath, '"',' ', fname];
%         system(systemCommand)
    end
end