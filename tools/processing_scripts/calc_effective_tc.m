function [approxTc, leftpeakloc, rightpeakloc, valleyloc] = calc_effective_tc(fname, draw)
%CALC_EFFECTIVE_TC Finds the double peak temperature in given dos file
%   Reads in a DoS File and will find the double peak location. Will also
%   draw out a set of plots around approxTc if draw is true

[params, g, e, en] = readDoSFile(fname);

tc = 1 /(log(1+sqrt(params.numStates)));
threshold = 0.000001;
approxTc = tc;
lower = tc;
if params.numStates > 5
  upper = tc+0.004; % Set to 0.005 for q=5, 0.004 for q>5
else
  upper = tc+0.005;
end
mid = (lower+upper)/2;

% Simple binary search for 40 iterations
for idx = 1:40
  % fprintf('Iteration: %d\n', idx);
  approxTc = mid;
  ft = g-e./mid;
  pe = exp(ft-max(ft));
  [peaks, locs] = findpeaks(pe, 'MinPeakDistance', 1);
  peaks = peaks(en(locs) > -1.8);
  locs = locs(en(locs) > -1.8);
  leftpeaks = peaks(en(locs) < -1.4);
  rightpeaks = peaks(en(locs) > -1.4);
  leftlocs = locs(en(locs) < -1.4);
  rightlocs = locs(en(locs) > -1.4);
  %[leftpeak, lidx] = max(peaks(en(locs) < -1.4));
  %[rightpeak, ridx] = max(peaks(en(locs) > -1.4));
  [leftpeak, lidx] = max(leftpeaks);
  [rightpeak, ridx] = max(rightpeaks);

   
  if isscalar(leftpeak) && isscalar(rightpeak)
    if (1-leftpeak < threshold) && (1- rightpeak < threshold)
      approxTc = mid;
      leftpeakloc = en(leftlocs(lidx));
      rightpeakloc = en(rightlocs(ridx));
      peupsidedown = -pe;
      [peaks1, locs1] = findpeaks(peupsidedown, 'MinPeakDistance', 1);
      mask = en(locs1) > leftpeakloc & en(locs1) < rightpeakloc;
      peaks1 = peaks1(mask);
      locs1 = locs1(mask);
      [valleypeak, vidx] = max(peaks1);
      valleyloc = en(locs(vidx));
      % fprintf('Found t, breaking\n');
      break;
    elseif leftpeak > rightpeak
      % fprintf('Increasing\n');
      lower = mid;
      mid = (lower + upper) / 2;
    else
      % fprintf('Decreasing\n');
      upper = mid;
      mid = (lower + upper) /2;
    end
  end
end

fprintf('T_c = %f, approxtc = %f\n', tc, approxTc);

if draw
  for ttest = approxTc-0.0004:0.0001:approxTc+0.0004
    %ttest=approxTc
    ft = g-e./ttest;
    pe = exp(ft-max(ft));
    figure;
    hold on;
    %titleStr = sprintf('P(T), q=%d, L=%d, Tc=%f, Approx Tc=%f\n', params.numStates, params.L, tc, approxTc);
    titleStr = sprintf('P(E,T), q=%d, L=%d, T=%f\n', params.numStates, params.L, ttest);
    title(titleStr);
    xlabel('-E/N');
    ylabel('P(E,T)');
    n = sprintf('T=%f', ttest);
    plot(en, pe, 'DisplayName', n);
    set(gca, 'YScale', 'log');
    ylim([1e-4, 1]);
    xlim([-2, 0]);
    hold off;
  end
end