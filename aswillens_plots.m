horz_vel = pivfilt.U;
vert_vel = -pivfilt.V;
%% Frame wise representation
% I = nanmean(sqrt(horz_vel.^2 + vert_vel.^2), 3);
% Ic = round( I/max(max(I))*64);
% Ic( Ic == 0) = 1;
% Ic(isnan(Ic)) = 1;
% colormap hsv;
% C = colormap;
% 
% log_PDI = mean(abs(BF_filtered),3); 
% log_PDI = imresize(log_PDI,[ kgrid.Ny , kgrid.Nx]);% should be done in beamformed dimensions
% 
 frames = size(horz_vel,3);
% f1 = floor(frames/3); f2 = floor(frames*2/3); f3 = frames - 1;
% 
% figure; 
% hax = [piv.X(1,1), piv.X(1,end)];
% vax = [piv.Z(1,1), piv.Z(1,end)];
% imagesc(0*log10(abs(log_PDI/max(log_PDI(:)))));
% colormap gray;
% hold on;
% 
% reso = abs(piv.X(1,:) - piv.X(1,2));
% lvquiv_mf = 1;
% lvquiv_sf = 2;
% 
%  hf1 = horz_vel(:,:,f1).*lvquiv_sf;
%  vf1 = vert_vel(:,:,f1).*lvquiv_sf;
% 
% 
% 
% for indsz = 1:lvquiv_mf:size(piv.X,1)
% for indsx = 1:lvquiv_mf:size(piv.X,2)
% lvquiv_quiv = quiver(piv.X(indsz,indsx),piv.Z(indsz,indsx),...
%     hf1(indsz,indsx), vf1(indsz,indsx),0,...
%     'Color',C(Ic(indsz,indsx),:),'maxheadsize',2);
% end
% end
% 
% hold off;
% xlabel('Lateral direction in [mm]')
% ylabel('Axial direction in [mm]')
% title(['Frame ' num2str(f1) ]);

%% Ground truth
% I = nanmean(sqrt(horz_vel.^2 + vert_vel.^2), 3);
% Ic = round( I/max(max(I))*64);
% Ic( Ic == 0) = 1;
% Ic(isnan(Ic)) = 1;
% colormap hsv;
% C = colormap;
% 
% figure; 
% hax = [scan.x(1), scan.x(end)];
% vax = [scan.z(1), scan.z(end)];
% imagesc(0*log10(abs(log_PDI/max(log_PDI(:)))));
% colormap gray;
% hold on;
% 
% lvquiv_mf = 1;
% lvquiv_sf = 1;
% 
% for indsz = 1:lvquiv_mf:size(piv.X,1);
% for indsx = 1:lvquiv_mf:size(piv.X,2);
% lvquiv_quiv = quiver(piv.X(indsz,indsx),piv.Z(indsz,indsx),...
%     0.*lvquiv_sf, 1.*lvquiv_sf,0,...
%     'Color',C(Ic(indsz,indsx),:),'maxheadsize',2);
% end
% end
% hold off;
% xlabel('Lateral direction in [mm]');
% ylabel('Axial direction in [mm]');
% title(['Frame ' num2str(f1) ]);

%% scatter plot

TF = ~isnan(vert_vel);
l_index = find(TF);
v_z_matrix = squeeze(vert_vel(~isnan(vert_vel)  ));
v_x_matrix = squeeze(horz_vel(~isnan(horz_vel)  ));
v_z = v_z_matrix(:);
v_x = v_x_matrix(:);

sv_z = -1.5; sv_x = 0;
sv_z = 0.15; sv_x = 0.2598;
sv_z = 0; sv_x = 0.3;
sv_z = 0.3; sv_x = 0.5196;

corrvalues = piv.corrMap;
corrvalues(corrvalues<0.7) = 0;
corrvalues_matrix = squeeze(corrvalues( l_index  ));
%corrvalues_array = (corrvalues_matrix(:)- min(corrvalues_matrix(:)))./(max(corrvalues_matrix(:))-min(corrvalues_matrix(:)));
corrvalues_array = corrvalues_matrix(:);
map = colormap('jet(10)');
% corr_norm = round(corrvalues_array.*64);
corr_norm = round(corrvalues_array.*10);
corr_norm(corr_norm == 0) = 1;
corr_map = zeros(length(corr_norm),3);
corr_map = map(corr_norm,:);

figure; 
plot(v_z, 'o');
hold on 
plot( sv_z.*ones(1, length(v_z)) );
hold off
legend('Measured axial velocity', 'True axial velocity');
xlabel('samples')
ylabel('Measured axial velocity [m/s]');

figure; 
plot(v_x, '*');
hold on 
plot( sv_x.*ones(1, length(v_z)) );
hold off
legend('Measured lateral velocity', 'True lateral velocity');
xlabel('samples')
ylabel('Measured lateral velocity [m/s]');

figure; 
scatter(v_x, v_z,[],  corr_map);
hold on
plot(sv_x, sv_z,'*');
line([0 sv_x], [0 sv_z]);
hold off
colormap jet;
colorbar ;
caxis([0.7 1]);
legend('measured velocity','ground velocity');
% hold on 
% plot( sv_x.*ones(1, length(v_z)) );
% hold off
% legend('Measured lateral velocity', 'True lateral velocity');
ylabel('Measured axial velocity [m/s]')
xlabel('Measured lateral velocity [m/s]');
xlim([-0.5 2]);
ylim([-0.3 0.5]); %0.08

%% mean absolute deviation for the frames
% remove values more than 3*std from the mean
% frames should be minimum 20 for outlier detection to work(one outlier)

m_h = nanmean(horz_vel, 3);
m_v = nanmean(vert_vel, 3);
std_h = nanstd(horz_vel,  3);
std_v = nanstd(vert_vel,  3);

shift_m_h = horz_vel - m_h;
shift_m_v = vert_vel - m_v;
out_horz_vel = horz_vel; out_vert_vel = vert_vel;

out_horz_vel(abs(shift_m_h) >= 2*std_h) = NaN;
out_vert_vel(abs(shift_m_v) >= 2*std_v) = NaN;

out_horz_vel = reshape(out_horz_vel, [], frames);
out_vert_vel = reshape(out_vert_vel, [], frames);
mad_horz = nanmean(abs(out_horz_vel - sv_x));
mad_vert = nanmean(abs(out_vert_vel - sv_z));

figure;
plot(mad_vert,'o-');
xlabel('frames');
ylabel('axial deviation');

figure;
plot(mad_horz,'o-');
xlabel('frames');
ylabel('lateral deviation');
%% Scatter plot with removal
out_v_z = squeeze(out_vert_vel(~isnan(out_vert_vel)  ));
out_v_x = squeeze(out_horz_vel(~isnan(out_horz_vel)  ));

figure; 
plot(out_v_z, 'o');
hold on 
plot( sv_z.*ones(1, length(out_v_z)) );
hold off
legend('Measured axial velocity', 'True axial velocity');
xlabel('samples')
ylabel('Measured axial velocity [m/s]');
title('Outliers removed');

figure; 
plot(out_v_x, '*');
hold on 
plot( sv_x.*ones(1, length(out_v_z)) );
hold off
legend('Measured lateral velocity', 'True lateral velocity');
xlabel('samples')
ylabel('Measured lateral velocity [m/s]');
title('Outliers removed');
%% Linear regression analysis
% mean, standand deviation, LSE, bias
% values more than 3*std from the mean are removed 

horz_vel_all = reshape(out_horz_vel, 1, []);
vert_vel_all = reshape(out_vert_vel, 1, []);

mean_all_h = nanmean(horz_vel_all);
mean_all_v = nanmean(vert_vel_all);

bias_h = mean_all_h - sv_x;
bias_v = mean_all_v - sv_z;

std_all_h = nanstd(horz_vel_all);
std_all_v = nanstd(vert_vel_all);

mse_h = mean_all_h.^2 + std_all_h.^2;
mse_v = mean_all_v.^2 + std_all_v.^2;

Values = [mean_all_v, bias_v, std_all_v, mse_v];

VarNames = { 'mean', 'bias', 'std','MSE'};
fprintf(1, '  %s\t\t%s\t\t\t%s\t\t%s\n', VarNames{:})
fprintf(1, '\t%06f\t%06f\t%06f\t%06f\n', Values');

Values = [mean_all_h, bias_h, std_all_h, mse_h];

VarNames = { 'mean', 'bias', 'std','MSE'};
fprintf(1, '  %s\t\t%s\t\t\t%s\t\t%s\n', VarNames{:})
fprintf(1, '\t%06f\t%06f\t%06f\t%06f\n', Values');