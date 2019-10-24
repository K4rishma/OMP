function [vector,hort,vert] = RobustOMP(image1, image2, sw, masked_indices) 
%image1(16x16)%image2(32x32) - already segmented in blocks
% keeping symmetric blocks
% image 1 has three diminsions
% remove all sorts of normalisation

%9/9 normalisation removed, sw =1 , normaisation should be changed
% checking on slanted vessel
b_size = size(image1, 1);
se_size = size(image2, 1);
kernels = size(image1,3);
 D = zeros(b_size.^2, b_size.^2, kernels);
% D = zeros(b_size.^2, 16, kernels);

% D_1 = karishma_D(image2, b_size);
D_1 = shifted_D(image2, b_size);
%D_1 = (D_1 - mean(D_1))./std(D_1); % this normalization has its importance in feature normalization in machine learning

% for k = 1: kernels
% D(:,:,k) = normc(D_1(:,:,k));
% end

D = D_1;

% implement for ill conditioned pseudoinverse
y = reshape(image1, b_size.^2 ,kernels );
%y = (y - mean(y))./std(y);

% Initialisation
itr = 1;
itr_max  = 1;
 x = zeros(b_size.^2, kernels, itr_max); 
% x = zeros(16, kernels, itr_max); 
rho = zeros(b_size.^2, kernels, itr_max); 
rho(:,:,itr) =  sw.*y;
disp = NaN(kernels,itr_max);

%ii = ~mask(squeeze(ss1(round(ksize(1)/2+1), round(ksize(2)/2+1), 1, :))); % image pixel indices of mask
kernel_indices = 1:1:kernels;


for itr = 1: 1: itr_max
    
    if itr == 1
        residue = sw.*y;
    else
        residue = rho(:,:,itr-1);
    end
    
%     for k = kernel_indices(1):kernel_indices(end)
for k = 1:kernels
        
        arg = abs(D(:,:,k)'*residue(:,k));
        %index = find(arg == nanmax(arg) & arg~=0.001.*ones(b_size.^2), 1); % maximum more than one ?
        index = find(arg == nanmax(arg) & nanmax(arg)>0,1);% arg~=0.001.*ones(14), 1);% non correct
        % check its 3d
        
        if (itr >= 2) & disp(k,itr-1) == index%& (find(disp(k,1:itr-1) == (index)))
            disp(k,itr) = index; 
            rho(:,k,itr)  = rho(:,k,itr-1); 
%             kernel_indices(k) = [];
            
        elseif ~isnan(index)
            disp(k,itr) = index;
            phi = D(:, disp(k,1:itr), k);
            %x(disp(k,1:itr),k,itr) = pinv(D(:, disp(k,1:itr), k))*y(:,k); % check here too
            x(disp(k,1:itr),k,itr) = (inv(phi.'*phi)*phi.')*y(:,k);
            weights = sw;% include local spatial tuning another function
            %rho(:,k,itr) = weights.*(y(:,k) - D(:, disp(k,1:itr), k)*x(disp(k,1:itr), k, itr)); %check if x is sparse
            
            rho(:,k,itr) = (eye(b_size.^2) -phi*inv(phi.'*phi)*phi.')*y(:,k);
        end
        
    end
    
end


%subpixel = reshape(x(disp,itr_max+1),  %check the indexing(velocity exp)

%mdisp = nanmean(disp,2);
[vert,hort] = ind2sub([b_size,b_size],disp);
x_value = zeros(kernels,itr_max);

for k = 1:kernels
    for itr = 1:itr_max
        
        x_value(k,itr) = x(disp(k,itr),k,itr_max)';
    end
end
 % normalize 
x_value = normr(x_value);
%change this and check with subpixel fit
%  Y_disp = nanmean((vert - b_size/2 - 1).*x_value,2);
% X_disp = nanmean((hort - b_size/2 -1).*x_value,2);

x_disp = zeros(kernels,1); y_disp = zeros(kernels,1);
for k = 1:kernels
cor = xcorr2(image2(:,:,k)./std2(image2(:,:,k)),image1(:,:,k)./std2(image1(:,:,k)));
% cor = xcorr2(image1(:,:,k),image2(:,:,k));
[y_disp(k),x_disp(k)] = find(cor == max(max(cor)) )   ;
end
vert = y_disp;
hort = x_disp;
% Y_disp = nanmean(-b_size + vert - 1 ,2);
% X_disp = nanmean((-b_size + hort - 1 ),2);
Y_disp = -b_size + vert ;
X_disp = -b_size + hort ;

vector = zeros(kernels,2);
vector = [X_disp(:),Y_disp(:)];
%vector(masked_indices,:) = NaN;

end