function tglT_m = display_dos_data(fnames, R, T, dosFile)

numr = length(R);
numt = length(T);

[params, g, e, en] = readDoSFile(dosFile);

tglT = zeros(numr, numt);
tglTgauss = zeros(numr, numt);
tglTb = zeros(numr, numt);
tglTs = zeros(numr, numt);
ipw = zeros(numr, numt);
Tvals = zeros(1, numt);
numE1=2*params.L^2+1-4
numE=2*32*32+1-4;
%numE=2*128*128+1-4;
%numE=2*64*64+1-4;
gteVals = zeros(numE, numt, numr);
gteBinaryVals = zeros(numE, numt, numr);
miVals = zeros(numE, numt, numr);
magnetisationVals = zeros(numE, numt, numr); % Before multiplication with P(E)
magCounts = zeros(numE, numt, numr);
nlogp = zeros(numE, numt, numr);
nlogpOrig = zeros(numE1, numt, numr);
ple = zeros(numE, numt, numr);
mag = zeros(numr, numt);
%mag_count = zeros(numE, numt, numr);
%[~, Eo, Ed, ~] = calc_effective_tc(dosFile, false);
[~, Eo, Ed, ~] = find_e(params.numStates);
[~, Eo, Ed, ~] = calc_effective_tc(dosFile, false);
tc= 1/log(1+sqrt(params.numStates));
for r = 1:numr
    for t = 1:numt
        draw = false;
        [paramsGTE, tglT(r,t), Tvals(:,t), gte, en, nlogp1, enOrig, nlogpOrig1, tglTgauss(r,t),tglTb(r,t), gteb, mag(r,t), magCount1, ple1, mag1, ipw(r,t), mi, tglTs(r,t)] = calc_dos_data(fnames{t, r}, dosFile, draw, Eo, Ed);
        
        gteVals(:,t, r)=gte;
        gteBinaryVals(:,t,r) = gteb;
        miVals(:,t,r)=mi;
        magnetisationVals(:,t,r) = mag1;
        magCounts(:,t,r) = magCount1;
        nlogp(:,t,r)= nlogp1;
        nlogpOrig(:,t,r)= nlogpOrig1;
        ple(:,t,r) = ple1;
    end
end
paramsGTE
ColorSet = varycolor(numt);

fprintf('Stats\n');
magCountSum = sum(magCounts,2);
magCountSum(magCountSum==0)=1; % /1 does nothing. Can only get 0 if magCounts are all 0 for that E
% value, thus 0/1 = 0 as desired
prop = bsxfun(@rdivide, magCounts, magCountSum);
% Need to calc the weighted average, not use prop
collapsedMag = sum(prop.*magnetisationVals, 2);
newMag = zeros(numr, numt);
newMagPLE = zeros(numr, numt);
for r = 1:numr
    for t =1:numt
        newMag(r,t) = sum(collapsedMag(:,1,r) .* nlogp(:, t, r)) / (sum(nlogp(:,t,r)));
        newMagPLE(r,t) = sum(collapsedMag(:,1,r) .* ple(:, t, r)) / (sum(ple(:,t,r)));
    end
end
figure;
hold on;
collapsedMag_m = mean(collapsedMag, 3);
plot(en, collapsedMag_m, '-r', 'DisplayName', 'collapsedMag');
title('Collapsed mag');
legendCell2 = cellstr(num2str(Tvals', 'T=%-f'));
ax2 = axes('YAxisLocation', 'Right');
hold on;
set(ax2, 'Color', 'none');
set(ax2, 'ColorOrder', ColorSet);
nlogp_m = mean(nlogp, 3);
plot(ax2, en, nlogp_m, 'DisplayName', legendCell2);
legend('show');
hold off;

% May need to fix the above to account for reps?
gteCountSum = sum(magCounts,2);
gteCountSum(gteCountSum==0)=1; % /1 does nothing. Can only get 0 if magCounts are all 0 for that E
% value, thus 0/1 = 0 as desired
prop = bsxfun(@rdivide, magCounts, gteCountSum);
% Need to calc the weighted average, not use prop
collapsedGte = sum(prop.*gteVals, 2);
% return
newGte = zeros(numr, numt);
for r = 1:numr
    for t =1:numt
        newGte(r,t) = sum(collapsedGte(:,1,r) .* nlogp(:, t, r)) / (sum(nlogp(:,t,r)));
    end
end
newGte = newGte .* log2(exp(1));
figure;
hold on;
collapsedGte_m = mean(collapsedGte, 3);
plot(en, collapsedGte_m, '-r', 'DisplayName', 'collapsedGTE');
title('Collapsed GTE');
legendCell2 = cellstr(num2str(Tvals', 'T=%-f'));
ax2 = axes('YAxisLocation', 'Right');
hold on;
set(ax2, 'Color', 'none');
set(ax2, 'ColorOrder', ColorSet);
size(en)
nlogp_m = mean(nlogp, 3);
plot(ax2, en, nlogp_m, 'DisplayName', legendCell2);
legend('show');
hold off;

tglT_m = mean(tglT, 1)';
tglT_s = std(tglT,0,1)'/sqrt(numr);
tglTgauss_m = mean(tglTgauss, 1)';
tglTgauss_s = std(tglTgauss,0,1)'/sqrt(numr);
tglTb_m = mean(tglTb, 1)';
tglTb_s = std(tglTb,0,1)'/sqrt(numr);
tglTs_m = mean(tglTs, 1)';
tglTs_s = std(tglTs,0,1)'/sqrt(numr);
tglT1_m = mean(newGte,1)';
tglT1_s = std(newGte,0,1)'/sqrt(numr);
ipw_m = mean(ipw, 1)';
ipw_s = std(ipw,0,1)'/sqrt(numr);

mag_m = mean(mag, 1)';
mag_s = std(mag,0,1)'/sqrt(numr);
mag1_m = mean(newMag, 1)';
mag1_s = std(newMag, 0, 1)'/sqrt(numr);
mag2_m = mean(newMagPLE, 1)';
mag2_s = std(newMagPLE, 0, 1)'/sqrt(numr);

figure;
hold on;
errorbar(Tvals, tglT_m, tglT_s, '-bo', 'DisplayName', 'GTE - DoS');
errorbar(Tvals, tglTgauss_m, tglTgauss_s, '-rs', 'DisplayName', 'GTE - Gaussian Approximation');
errorbar(Tvals, tglTb_m, tglTb_s, '-c>', 'DisplayName', 'GTE - Binary');
errorbar(Tvals, tglTs_m, tglTs_s, '-c>', 'DisplayName', 'GTE - Energy - Sweep');
errorbar(Tvals, tglT1_m, tglT1_s, '--gx', 'DisplayName', 'GTE - DoS - Collapsed');
errorbar(Tvals, ipw_m, ipw_s, '-md', 'DisplayName', 'MI - DoS (sweep)');
line([tc tc], ylim, 'Color', 'black', 'LineStyle','--');
titleStr = sprintf('DoS Metrics: q=%d, L=%d, U=%d, S=0, T_c=%f', paramsGTE.numStates, paramsGTE.L, paramsGTE.U, tc);
title(titleStr);
xlabel('T');
ylabel('Bits');
legend('show');

ax2 = axes('YAxisLocation', 'Right');
hold on;
set(ax2, 'Color', 'none');
errorbar(Tvals, mag_m, mag_s, '-gd', 'DisplayName', 'Magnetisation');
errorbar(Tvals, mag1_m, mag1_s, '--bs', 'DisplayName', 'Magnetisation (Collapsed)');
%errorbar(Tvals, mag2_m, mag2_s, '--bs', 'DisplayName', 'Magnetisation (Collapsed, PLE)');
ylabel('Magnetisation');
hold off;

figure;

set(gca, 'ColorOrder', ColorSet);
hold on;
gteVals_m = mean(gteVals,3);
plot(en, gteVals_m); % All
titleStr = sprintf('GTE before multiplying by P(E): q=%d, L=%d, U=%d, S=0, T_c=%f', paramsGTE.numStates, paramsGTE.L, paramsGTE.U, tc);
title(titleStr);
legendCell = cellstr(num2str(Tvals', 'T=%-f'));
legendCell1 = cellstr(num2str(Tvals', 'T=%-f (Gauss)'));
legend(legendCell);

ax2 = axes('YAxisLocation', 'Right');
hold on;
set(ax2, 'Color', 'none');
set(ax2, 'ColorOrder', ColorSet);
nlogp_m = mean(nlogp, 3);
nlogpOrig_m = mean(nlogpOrig, 3);
ple_m = mean(ple, 3);
plot(ax2, en, nlogp_m, 'DisplayName', legendCell);
plot(ax2, en, ple_m, 'DisplayName', legendCell1);
plot(ax2, enOrig, nlogpOrig_m, 'DisplayName', legendCell);


figure;
ColorSet = varycolor(numt);
set(gca, 'ColorOrder', ColorSet);
hold on;
titleStr = sprintf('Magnetisation before multiplying by P(E): q=%d, L=%d, U=%d, S=0, T_c=%f', paramsGTE.numStates, paramsGTE.L, paramsGTE.U, tc);
title(titleStr);
magnetisationVals_m = mean(magnetisationVals,3);

plot(en, magnetisationVals_m); % All
legendCell = cellstr(num2str(Tvals', 'T=%-f'));
legend(legendCell);

ax2 = axes('YAxisLocation', 'Right');
hold on;
set(ax2, 'Color', 'none');
set(ax2, 'ColorOrder', ColorSet);
plot(ax2, en, nlogp_m, 'DisplayName', legendCell);
hold off;
legend('show');


figure;
set(gca, 'ColorOrder', ColorSet);
hold on;
gteBinaryVals_m = mean(gteBinaryVals,3);
plot(en, gteBinaryVals_m); % All
titleStr = sprintf('GTE Binary before multiplying by P(E): q=%d, L=%d, U=%d, S=0, T_c=%f', paramsGTE.numStates, paramsGTE.L, paramsGTE.U, tc);
title(titleStr);
legendCell = cellstr(num2str(Tvals', 'T=%-f'));
legend(legendCell);

ax2 = axes('YAxisLocation', 'Right');
hold on;
set(ax2, 'Color', 'none');
set(ax2, 'ColorOrder', ColorSet);
plot(ax2, en, nlogp_m, 'DisplayName', legendCell);


figure;
set(gca, 'ColorOrder', ColorSet);
hold on;
miVals_m = mean(miVals,3);
plot(en, miVals_m); % All
titleStr = sprintf('MI before multiplying by P(E): q=%d, L=%d, U=%d, S=0, T_c=%f', paramsGTE.numStates, paramsGTE.L, paramsGTE.U, tc);
title(titleStr);
legendCell = cellstr(num2str(Tvals', 'T=%-f'));
legend(legendCell);

ax2 = axes('YAxisLocation', 'Right');
hold on;
set(ax2, 'Color', 'none');
set(ax2, 'ColorOrder', ColorSet);
plot(ax2, en, nlogp_m, 'DisplayName', legendCell);