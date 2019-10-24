function [tprof_us, tstd_us, tax] = PlotProfile( us_x, us_y, us_u, us_v, tRes, mode )
addpath(genpath(getenv('PIVPATH')));
if nargin < 6
    mode = 'mag';  % can also choose u, v
end
%% Comparison App
us_mag = sqrt(us_u.^2 + us_v.^2);
us_mip = squeeze(nanmax(us_mag,3));
usvec0 = [];
tax = [1:size(us_u,3)].*tRes;
fig = figure();
magmin = min(us_mip(:));
magmax = max(us_mip(:));
% max mip image
ax1 = subplot(2,2,1);
magim = imagesc(us_x(1,:), us_y(:,1), us_mip, [magmin, magmax]);
axis image
axis ij
title('ePIV')
ylabel('Depth [mm]')
xlabel('Azimuth [mm')
ust = text(max(us_x(:))+10,min(us_y(:)),num2str(norm(diff(usvec0,1,1))));
while ishandle(fig)
    ust.String = num2str(norm(diff(usvec0,1,1)));  
    hus = imline(ax1,usvec0);
    usvec = wait(hus);
    hus.delete;
    usvecs = sort(usvec);
    if ~isequal(usvec,usvec0) && ishandle(fig)
        yr = mean(diff(us_y(:,1)));
        xr = mean(diff(us_x(1,:)));
        us_indsx = find((us_x(1,:)>=(usvecs(1,1)-xr)) & (us_x(1,:)<=(usvecs(2,1)+xr)));
        us_indsy = find((us_y(:,1)>=(usvecs(1,2)-yr)) & (us_y(:,1)<=(usvecs(2,2)+yr)));
        us_u_ol = us_u(us_indsy,us_indsx,:);
        us_v_ol = us_v(us_indsy,us_indsx,:);
        xsec_us_u = squeeze(nanmean(us_u_ol,1));
        xsec_us_v = squeeze(nanmean(us_v_ol,1));
        xsec_us_mag = sqrt(xsec_us_u.^2 + xsec_us_v.^2);
        switch mode
            case 'mag'
                tprof_us = squeeze(nanmean(xsec_us_mag,1));
                tstd_us = squeeze(nanstd(xsec_us_mag,1));
            case 'u'
                tprof_us = squeeze(nanmean(xsec_us_u,1));
                tstd_us = squeeze(nanstd(xsec_us_u,1));
            case 'v'
                tprof_us = squeeze(nanmean(xsec_us_v,1));
                tstd_us = squeeze(nanstd(xsec_us_v,1));
            case 'mag_signed'
                tprof_us = squeeze(nanmean(xsec_us_mag,1));
                tprof_us = tprof_us.*sign(nanmean(xsec_us_v,1)).*-1;
                tstd_us = squeeze(nanstd(xsec_us_mag,1));
            case 'norm'
                lnormal = flip(-1.*diff(usvec,1,1));
                usvecs = zeros([2, size(xsec_us_u)]);
                usvecs(1,:,:) = xsec_us_u;
                usvecs(2,:,:) = xsec_us_v;
                usvecs_proj = dot(usvecs,repmat(lnormal',1, size(usvecs,2), size(usvecs,3)),1)./(norm(lnormal));
                tprof_us = squeeze(nanmean(usvecs_proj,2));
                tstd_us = squeeze(nanstd(usvecs_proj,2));  
        end
        usvec0 = usvec;
    end
    if ishandle(fig)
        ust.String = num2str(norm(diff(usvec0,1,1)));
        %plot comparison
        subplot(2,2,[3,4])
        shadedErrorBar(tax, tprof_us, tstd_us, {'Color',[0.8500, 0.3250, 0.0980]}, 1)
        title('Velocity time profile comparison')
        xlabel('Time [s]')
        ylabel('Velocity [m/s]')
        pause(0.1)  
    end
end


end

