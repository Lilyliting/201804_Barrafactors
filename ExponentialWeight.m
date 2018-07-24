function w = ExponentialWeight(T, hl)
t = (0:T)';
w = flipud(exp(-log(2)/hl*t));
end