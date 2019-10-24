function [flow_grid] = get_SAVFI_grid(filepath)
    %%% determine which SAVFI phantom is being used based on file
    %%% structure
    savfi_eval = 'SAVFI/reference_grids';
    phantoms = dir(savfi_eval);
    phantoms = phantoms(3:end);
    % check what type of phantom is being used (assumes filepath is
    % preseved more or less from challenge)
    for i = 1:length(phantoms)  % check phantom type from filepath
        if contains(filepath, phantoms(i).name)
            grid_path=fullfile(savfi_eval, phantoms(i).name);
            % first check if cfd model - has training and testing
            % subfolder
            if contains(filepath,'carotid_bifurcation')
                datatype = dir(grid_path);
                datatype = datatype(3:end);
                for j = 1:length(datatype)
                    if contains(filepath, datatype(j).name)
                        grid_path=fullfile(grid_path, datatype(j).name);
                        break
                    end
                end
            end
            % then check whether phantom is simulation or measurement   
            phantom_type = dir(grid_path);
            phantom_type = phantom_type(3:end);
            for j = 1:length(phantom_type)
                if contains(filepath, phantom_type(j).name)
                    grid_path=fullfile(grid_path, phantom_type(j).name);
                    break
                end
            end
            break
        end
    end
    A = load(fullfile(grid_path, 'ref_grid.mat'));
    flow_grid = A.flow_grid;
end