function [F,Err] = MultiFactorsRegression(X, Y, W, mode)
% 


if size(Y,1) ~= 1
    Y = Y';
    X = X';
end
if size(W,1) ~= 1
    W = W';
end

%N = size(Y,2);
f = size(X,1);
A = zeros(f,f);
B = zeros(f,1);
for ii = 1 : f
    B(ii,1) = nansum(W .* Y .* X(ii,:));
    for jj = ii : f
        A(ii,jj) = nansum(W .* X(ii,:) .* X(jj,:));
        if ii ~= jj
            A(jj,ii) = nansum(W .* X(ii,:) .* X(jj,:));
        end
    end
end
if strcmp(mode, 'real') == 1
    F = real((A \ B)');
else
    F = (A \ B)';
end
Err = Y - F * X;

end