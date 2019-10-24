% take mean of images 
%% angular position of the pixels in the block/window
% Input ksize
Y = repmat([ksize(1): -1: -(ksize(1) - 1)]', 1, ksize(2));
X = repmat([-ksize(2): 1: (ksize(2) - 1)], ksize(1), 1);

angular_pos = atan2d(Y,X);
%%% Weighing finctions %%%
s = 1;









