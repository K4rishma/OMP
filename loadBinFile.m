function [Frames] = loadBinFile(infopath, binpath)
%loadBinFile loads complex IQ data from binary file.
try
    tic; fprintf('Loading Data\n');
    load(infopath); %load info file;
    fid = fopen(binpath, 'r');
    Frames = fread(fid,prod(SaveShape),['*' DataType]);
    fclose(fid);
    Frames = Frames(1:SaveShape(1))+1i*Frames(1+SaveShape(1):prod(SaveShape));
    Frames = reshape(Frames,DataShape);
    t = toc; fprintf('loading file took %f seconds \n',t);
catch
    fprintf('Could not find file: %s, \n',binpath);
    return
end

end

