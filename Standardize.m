function [output] = Standardize(mtx)
% 去极值并标准化
output = mtx;
tmpmean = nanmean(mtx,2); 
tmpstd = nanstd(mtx,1,2);
maxline = tmpmean+3*tmpstd;
minline = tmpmean-3*tmpstd;
for ii = 1:size(mtx,1)
%nanpos = 1-isnan(mtx(ii,:));
mtxii = mtx(ii,:);
mtxii(mtxii>maxline(ii)) = maxline(ii);
mtxii(mtxii<minline(ii)) = minline(ii);
mtxii = (mtxii-nanmean(mtxii)) ./ nanstd(mtxii);
output(ii,:) = mtxii;
end

return
end