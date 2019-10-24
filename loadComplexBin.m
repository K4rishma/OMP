function [Frames] = loadComplexBin(binpath, SaveShape, DataType, DataShape)
    try
            tic; fprintf('Loading Data\n');
            fid = fopen(binpath, 'r');
            Frames = fread(fid,prod(SaveShape),['*' DataType]);
            fclose(fid);
            Frames = Frames(1:SaveShape(1))+1i*Frames(1+SaveShape(1):prod(SaveShape));
            Frames = reshape(Frames,DataShape);
            t = toc; fprintf('loading file took %f seconds \n',t);
        catch
            fprintf('Could not find file: %s, \n',FILENAME);
            return
    end 
end