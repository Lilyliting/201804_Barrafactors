function [output] = nandatafill(mtx)
A = mtx;
%A=stkGrossProfitMarginMtrx;
% nanpos1 = isnan(A(1,:));
% A(1,nanpos1) = 0;

for ii = 2:size(A,1)
    zeropos = A(ii,:)==0;
    nanpos = isnan(A(ii,:));
    A(ii,zeropos) = A(ii-1,zeropos);
    A(ii,nanpos) = A(ii-1,nanpos);
end
output = A;
return;
end