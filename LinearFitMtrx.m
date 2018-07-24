function [alpha,beta] = LinearFitMtrx(xt,yt,wExp)
% xt for x, yt for y, wExp for ExponentialWeight
wExpExpand = repmat(wExp,1,size(xt,2));
a11 = repmat(sum(wExp),1,size(xt,2));
a12 = sum(wExpExpand .* xt,1);
a21 = a12;
a22 = sum(wExpExpand .* (xt.^2),1);
b1 = sum(wExpExpand .* yt,1);
b2 = sum(wExpExpand .* xt .* yt,1);

alpha = (b1 .* a22 - b2 .* a12) ./ (a11 .* a22 - a21 .* a12);
beta = (b1 .* a21 - b2 .* a11) ./ (a12 .* a21 - a11 .* a22); 

end