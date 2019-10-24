function weighd_corrMaps = directional_constraint(ss1,corrMaps,angles)
zero_angles = zeros(size(angles)+1);
zero_angles(1:end-1,1:end-1) =  angles; % padding with zero
clear angles
angles = zero_angles;
angles = medfilt2(angles,[10 10]);
%  figure; imagesc(abs(angles))
ss1_angle = ss1(:,:,1,:);
angles_seg = angles(ss1_angle);
kernels = size(ss1_angle,4);
angles_seg_vect = reshape(angles_seg,[],kernels);
med_angles = nanmedian(angles_seg_vect);% nanmeadian nanmean

% re = reshape(med_angles, 49, 49);
% figure; imagesc(re);

dw = zeros(size(corrMaps));
dim1 = size(corrMaps,1);
dim2 = size(corrMaps,2);
temp = zeros(dim1,dim2);
temp(dim1/2-4:dim1/2+4,dim1/2-4:dim1/2+4) = 1;
temp(dim1/2:dim1/2+1,:) = 1;
base = temp;

for theta = -10:4:10
    temp_rotate = imrotate(temp,theta,'nearest','crop');
    base = base + temp_rotate;
end


base = base > 0; %binary scaling

% base(base==0) = [];
% debug = NaN(16,16);
% debug = NaN; % taking the affect of negative values(not correct)
% debug(9:16,1:9) = 1;
for nk = 1:kernels
    theta = med_angles(nk);
%     theta = nanmedian(nanmedian(angles_seg(:,:,nk)));
    if isnan(theta)
        dw(:,:,nk) = 0;
    else
        weights = imrotate(base,theta,'nearest','crop');
        dw(:,:,nk) = weights;
    end
%     if theta < 10 && theta > 5
%      
%     end
% dw(:,:,nk) = base;    
end
% mod  = corrMaps(:,:,nk)
% mod  = corrMaps(:,:,nk);[r,c]=max(mod(:))

% [r,c]=max(mod(:).*base(:))
weighd_corrMaps = corrMaps.*dw;
end