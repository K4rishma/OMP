function [s_svd_all,U,Sigma,V] = svdfilter(bf,min_reg_val,max_reg_val)

% This function computes the Doppler signal of the beamformed data by
% singular value decomposition. it follows the publication doi:10.1109/TMI.2015.2428634

% rejection_value represents the cutoff for the singular vectors. e.g. if
% set to 10, the first 9 singular vectors are cut out and the signal is
% only computed with the decomposed data from 10 to the ensemble maximum number 

% bf should have dimension: samples // lines // ensembles // angles (optional)

% Jason Voorneveld 10/02/2016

for na = 1:size(bf,4)
    bf_temp = squeeze(bf(:,:,:,na));
    S = reshape(bf_temp,[size(bf_temp,1)*size(bf_temp,2) size(bf_temp,3)]);
    S = S-mean(S,2);
    [U,Sigma,V] = svdecon(S);
    remove_arr = [1:min_reg_val-1,size(Sigma,1)-max_reg_val:size(Sigma,1)];
    Sigma(remove_arr,remove_arr) = 0;
    s_svd = U*Sigma*V';
    s_svd_all(:,:,:,na) = reshape(s_svd,size(bf_temp));    
end

end
               

