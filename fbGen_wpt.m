%% Stationary wavelet packet warped filterbank generation
J = 4;
pow = 2^J;
% wavelet name
wavelet = 'db10';
% wavelet = 'dmey';
% wavelet filters
[lo,hi,lo_r,hi_r] = wfilters(wavelet);
% lscdf = liftwave('cdf5.5');
% [lo,hi,lo_r,hi_r] = ls2filt(lscdf);
centrePt = 4096;
kDelta = zeros(8193, 1);
kDelta(centrePt) = 1;
kDelta = [zeros(pow * 4, 1); kDelta; zeros(pow * 4, 1)];
[~, pos] = max(kDelta);
actualTransformationLen = length(kDelta);
needPadding1 = 0;
if rem(actualTransformationLen,pow)>0
    needPadding1 = ceil(actualTransformationLen/pow) * pow - actualTransformationLen;
end
signal = [kDelta; zeros(needPadding1, 1)];
%% Forward SWPT
wpt = swpt(signal, lo, hi, J);
firstIdx = zeros(pow, 1);
for ch = 1 : pow
    for idx = 1 : size(wpt, 2)
        if abs(wpt(ch, idx)) > 0
            firstIdx(ch) = idx;
            break;
        end
    end
end
firstIdxFin = min(firstIdx);
lastIdx = zeros(pow, 1);
for ch = 1 : pow
    for idx = firstIdxFin : size(wpt, 2)
        if abs(wpt(ch, idx)) > 0
            lastIdx(ch) = idx;
        end
    end
end
lastIdxFin = max(lastIdx);
analysis_filters = wpt(:, firstIdxFin : lastIdxFin);
wpt2 = swpt(signal, lo_r, hi_r, J) / pow;
firstIdx2 = zeros(pow, 1);
for ch = 1 : pow
    for idx = 1 : size(wpt2, 2)
        if abs(wpt2(ch, idx)) > 0
            firstIdx2(ch) = idx;
            break;
        end
    end
end
firstIdx2Fin = min(firstIdx2);
lastIdx2 = zeros(pow, 1);
for ch = 1 : pow
    for idx = firstIdx2Fin : size(wpt2, 2)
        if abs(wpt2(ch, idx)) > 0
            lastIdx2(ch) = idx;
        end
    end
end
lastIdx2Fin = max(lastIdx2);
synthesis_filters = wpt2(:, firstIdxFin : lastIdxFin);
%% Warping
zpd = 2048;
resAna = zeros(size(analysis_filters, 1), size(analysis_filters, 2) + zpd + 1);
resSyn = zeros(size(synthesis_filters, 1), size(synthesis_filters, 2) + zpd + 1);
for idx = 1 : size(analysis_filters, 1)
    padded = [analysis_filters(idx, :), zeros(1, zpd)];
    resAna(idx, :) = warp_impres(padded, -0.4, size(padded, 2));
    padded = [synthesis_filters(idx, :), zeros(1, zpd)];
    resSyn(idx, :) = warp_impres(padded, -0.4, size(padded, 2));
end
%% Plot filterbank
[m,n]=size(resAna);
figure;
for i=1:m
    [H(i,:),W]=freqz(resAna(i,:),1);
    hold on
    plot(W,20*log10(abs(H(i,:))),'linewidth',2)
    axis([0 pi -100 18])
end


xlabel('Normalized Frequency','fontsize',12)
ylabel('Amplitude Response','fontsize',12)
hold off
%% Perfect reconstruction test
ch1 = conv(analysis_filters(1, :), synthesis_filters(1, :));
ch1W = conv(resAna(1, :), resSyn(1, :));
for idx = 2 : size(analysis_filters, 1)
    ch1 = ch1 + conv(analysis_filters(idx, :), synthesis_filters(idx, :));
    ch1W = ch1W + conv(resAna(idx, :), resSyn(idx, :));
end
fvtool(ch1);
fvtool(ch1W);