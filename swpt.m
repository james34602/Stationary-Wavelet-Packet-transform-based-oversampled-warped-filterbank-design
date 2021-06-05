function [wpt, pseudoFreq] = swpt(signal, lo, hi, J)
% Create array to hold wavelet packets and packet levels
% Initially create full tree
cfs = zeros(2 ^ (J + 1) - 2, length(signal));
cfs(1,:) = signal;
packetlevels = zeros(sum(2 .^ (1 : J)), 1);
times = 2 .^ (1 : J);
for ii = 1 : J
    for jj = 1 : times(ii)
        packetlevels(times(ii) - times(1) + jj) = ii;
    end
end
% Indices for first level
Idx = 1:2;
for kk = 1 : J
    index = 0;
    %Determine first packet for a given level
    jj = (2 ^ kk) - 1;
    if (kk > 1)
        Idx = find(packetlevels == (kk - 1));
    end
    for nn = 0 : ((2 ^ kk) / 2 - 1)
        index = index + 1;
        X = cfs(Idx(index), :);
        [appoxi, detl] = swtSingleLevel(X, lo, hi, kk - 1);
        if mod(nn, 2) == 0 %iseven(nn)
            cfs(jj + 2 * nn + 1, :) = detl;
            cfs(jj + 2 * nn, :) = appoxi;
        else
            cfs(jj + 2 * nn, :) = detl;
            cfs(jj + 2 * nn + 1, :) = appoxi;
        end
    end
end
wpt = cfs(end-2^J+1:end,:);
% Generating vector of frequencies by level if needed
packetlevels = packetlevels(end-2^J+1:end);
df = 1./2.^(packetlevels+1);
df(1) = df(1) / 2;
pseudoFreq = cumsum(df);
end
function [appoxi, detl] = swtSingleLevel(x,lo,hi,J)
% Use row vector.
x = x(:)';
s = length(x);
% Compute stationary wavelet coefficients.
detl = zeros(1,s);
appoxi = zeros(1,s);
temp_lo = lo;
temp_hi = hi;
for ii = 1:J
    % upsample filters.
    temp_lo = dyadup(temp_lo,0,1);
    temp_hi = dyadup(temp_hi,0,1);
end
% Extension.
lf = length(temp_lo);
xExt  = [zeros(1, lf/2), x, zeros(1, lf/2)];
% Decomposition.
hRes = conv(xExt,temp_hi);
lRes = conv(xExt,temp_lo);
detl(:) = hRes(lf + 1 : s + lf);
appoxi(:) = lRes(lf + 1 : s + lf);
end