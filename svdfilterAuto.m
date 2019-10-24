function [s_svd,U, Sigma, V, lambdas1, lambdas2] = svdfilterAuto(IQ,Tthresh,rem_noise,Nthresh)

% This function computes the Doppler signal of the beamformed data by
% singular value decomposition. it follows the publication doi:10.1109/TMI.2015.2428634

% Automatic cutoff values based on eigen energy for lower bound and mean
% doppler frequency for upper bound

% bf should have dimension: axial//azimuth//ensembles//(optional)angles

% input variables:
% IQ        - input image set - dimensions = [axial, azimuth,ensembles,(ang)] 
%             (z, x, t, angles) - angles are optional
% Tthresh   - the slope threshold used for the automatic low rank cutoff
%             detection - a value of between 0.9 to 0.99 works (i found 
%             0.99 worked the best. 
% rem_noise - boolean: 1 = do high rank truncation to remove noise. 0 = do
%             not
% Nthresh   - Value between 0 and 1. Percentage of Nyquist slow time
%             sampling rate (frame rate) to set high rank cutoff.
% Jason Voorneveld 22/02/2017

if nargin < 2
    Tthresh = 0.95;
    rem_noise=1;
    Nthresh = 0.9;
end
if nargin < 3
    rem_noise=1;
    Nthresh = 0.9;
end
if nargin < 4
    Nthresh = 0.9;
end


s_svd = zeros(size(IQ));  % init svd filtered frames
lambdas1 = zeros(size(IQ,4),1); % each angle has different lambda
lambdas2 = zeros(size(IQ,4),1);
for ang = 1:size(IQ,4)  %loop through angles
    bf = squeeze(IQ(:,:,:,ang)); % work on one angle at a time
    S = reshape(bf,[size(bf,1)*size(bf,2) size(bf,3)]); % casorati     
    [U,Sigma,V] = svdecon(S); % fast SVD
    % now determine ratio of successive singular values
    dRatio = diag(Sigma(2:end,2:end))./diag(Sigma(1:end-1,1:end-1));
    dRatio = movmean(dRatio, ceil(size(bf,3)/10)); % smooth a bit
    lambda1 = find(dRatio > Tthresh,1);
    if rem_noise  % if high rank truncation too
        for vi = size(V,2):-1:1 % compute doppler frequency using autocorr method
            auto1 = sum(V(1:end-1,vi).*conj(V(2:end,vi)))/(size(V,1)-1);
            dopV(vi) = angle(auto1)/pi;
        end
        dopV = moving(abs(dopV),3,'median');  %smooth a bit
        lambda2 = find(dopV > Nthresh,1);
    else
        lambda2 = size(Sigma,1)-1; %otherwise keep all non-low rank
    end
    remove_arr = [1:lambda1-1,lambda2+1:size(Sigma,1)]; 
    Sigma(remove_arr,remove_arr) = 0; % truncate sigma
    s_svd(:,:,:,ang) = reshape(U*Sigma*V',size(bf)); % recompute
    lambdas1(ang) = lambda1; % save values 
    lambdas2(ang) = lambda2;
end
end
               

