function [ u, v ] = std_outlier( u, v, thresh )
%std_outlier Replaces outliers with NaNs using a standard deviation threshold

meanu=nanmean(nanmean(u));
meanv=nanmean(nanmean(v));
std2u=nanstd(reshape(u,size(u,1)*size(u,2),1));
std2v=nanstd(reshape(v,size(v,1)*size(v,2),1));
minvalu=meanu-thresh*std2u;
maxvalu=meanu+thresh*std2u;
minvalv=meanv-thresh*std2v;
maxvalv=meanv+thresh*std2v;
u(u<minvalu)=NaN;
u(u>maxvalu)=NaN;
v(v<minvalv)=NaN;
v(v>maxvalv)=NaN;

end

