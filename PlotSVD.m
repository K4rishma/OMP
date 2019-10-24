function [sfig, vfig] = PlotSVD(U, Sigma,V, PRF, c, f0)
%c = speed of sound
%f0 = center frequency of pulse
%PRF = pulse repetition frequency
if nargin < 4
    scale = false;
elseif nargin < 6
    scale = true;
    velscale = false;
end
    %% Plot Sigma
    S = max(Sigma);
    sfig = figure();
    plot(S)
    title('Eigenvalue Energy')
    ylabel('Energy')
    xlabel('Singular Value')
    %% Plot PSD of V
    if scale
        [psdV, fax] = periodogram(V, [], [], PRF);
        if velscale
            yax = fax.*c/(2*f0);
            ystr = 'Mode Velocity in m/s';
        else
            yax = fax;
            ystr = 'Mode frequency [Hz]'
        end
    else
        [psdV, fax] = periodogram(V);
        yax = fax/pi;
        ystr = 'Proportion of Nyquist Frequency';
    end
    psdV = 20*log10(psdV);
    psdV = psdV - max(psdV(:));
    xax = [1, size(V,2)];
    vfig = figure();
    imagesc(xax,yax,psdV,[-120 0])
    colormap('parula')
    title('PSD Map')
    xlabel('Singular Value')
    ylabel(ystr)
    cb = colorbar;
    ylabel(cb, 'PSD [dB]');
    
    %% Plot Spatial Similarity Matrix
    C = cov2corr(cov(abs(U)));
    C = C - diag(diag(C));
    ussfig = figure();
    imagesc(C)
    colormap('parula')
    colorbar()
    cb = colorbar;
    ylabel(cb, 'Correlation');
    title('Spatial Similarity Matrix')
end