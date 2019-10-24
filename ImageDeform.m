function [image_deformed, utable, vtable] = ImageDeform(image, xtable, ytable, utable, ...
    vtable,n_el, mini, maxi, step, ksize)
%ImageDeform Deforms the second image when using iterative PIV with window
%deformation.
%INPUTS:
%   image - Set of 2D images (MxNxEnsembles)
%   xtable, ytable - PIV grid
%   utable, vtable - Previous iterations PIV results
%   n_el, mini, maxi, step, ksize - this iterations grid parameters from
%   get_extents function.
%OUTPUTS:
% image_deformed - Set of 2D images (MxNxEnsembles) deformed by previous
% iterations displacement grid.
% utable - new sized utable (interpolated)
% vtable - new sized vtable (interpolated)
    % make copy of old grid. we need it later
    xtable_old=xtable;  
    ytable_old=ytable;
    % make new grid (might be a different size)
    xtable = repmat((mini(2):step(2):maxi(2)), n_el(1), 1) + ksize(2)/2;
    ytable = repmat((mini(1):step(1):maxi(1))', 1, n_el(2)) + ksize(1)/2;
    % interpolate displacements onto new grid
    utable=interp2(xtable_old,ytable_old,utable,xtable,ytable,'LINEAR');
    vtable=interp2(xtable_old,ytable_old,vtable,xtable,ytable,'LINEAR');
    % find number of extrapolations required - this is why only multiples
    % of 2 are allowed for step size calculation
    n_extrap = double(ceil(ksize./step)/2); 
    % now pad array of new utable (replicate edge values)
    U= padarray(utable, n_extrap, 'replicate');
    V= padarray(vtable, n_extrap, 'replicate');  
    %add 1 line around image for border regions... linear extrap  
    firstlinex=xtable(1,:);
    firstlinex_intp=interp1(1:1:size(firstlinex,2),firstlinex,1-n_extrap(2):1:size(firstlinex,2)+n_extrap(2),'linear','extrap');
    X=repmat(firstlinex_intp,size(xtable,1)+n_extrap(2)*2,1);
    firstliney=ytable(:,1);
    firstliney_intp=interp1(1:1:size(firstliney,1),firstliney,1-n_extrap(1):1:size(firstliney,1)+n_extrap(1),'linear','extrap')';
    Y=repmat(firstliney_intp,1,size(ytable,2)+n_extrap(1)*2);
    % now create full grid of all pixels in the image
    [X1, Y1] = meshgrid(X(1,1):1:X(1,end)-1, Y(1,1):1:Y(end,1)-1);
    % interpolate displacements at each pixel of the image
    U1 = interp2(X,Y,U,X1,Y1,'LINEAR');
    V1 = interp2(X,Y,V,X1,Y1,'LINEAR');
    image_deformed = zeros(size(X1), class(image)); % init image matrix
    xinds = 1:size(image,2);  % 
    yinds = (1:size(image,1))';
    for i = 1:size(image,3)
        image_deformed(:,:,i) = interp2(xinds,yinds,image(:,:,i),...
            X1+U1,Y1+V1,Opts.imDeform); %linear is 3x faster and looks ok...
    end

end

