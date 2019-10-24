function [sw,med_angles] =  spatialDistribution(angles, Weight_fnc, angular_pos,ss1)
%% To Do
% corrent for angles = 0
%%
%  figure; imagesc(abs(angles))
ss1_angle = ss1(:,:,1,:);
angles_seg = angles(ss1_angle);
kernels = size(ss1_angle,4);
angles_seg_vect = reshape(angles_seg,[],kernels);
med_angles = nanmedian(angles_seg_vect);% nanmeadian nanmean


height = size(angular_pos,1);
width = size(angular_pos,2);
sw = zeros(height, width, kernels);

for k = 1 : kernels
angle_shift = med_angles(k);
% if ~isnan(angle_shift)
%     pause()
% end
for i = 1:height
    
    for j = 1:width
        
        if angular_pos(i,j)>90 && angular_pos(i,j)<=180
            
            theta = (angular_pos(i,j)-180 - angle_shift);
            
            if theta < -90 
            theta = theta + 180 ;
            end
            sw(i,j,k) = Weight_fnc(theta);
        
        elseif angular_pos(i,j)>-180 && angular_pos(i,j)<=-90
            
            theta = (angular_pos(i,j)+180 - angle_shift);
            if theta >= 90 
            theta = theta - 180 ;
            end
            sw(i,j,k) = Weight_fnc(theta);
            
        
        else
            
              theta = (angular_pos(i,j) - angle_shift);
            
            if theta < -90 
            theta = theta + 180 ;
            end
            
            if theta >= 90 
            theta = theta - 180 ;
            end
            
            sw(i,j,k) = Weight_fnc(theta);
        
        end
        
    end
    
end

sw = min(sw,1);
sw = max(sw,0);
 sw((end-1)/2+1,(end-1)/2+1,:) = 1;
sw(sw==Inf) = 1;
%sw = reshape(sw,[],kernels);
% figure; imagesc(sw);colorbar;



end