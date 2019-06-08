function [tvals, gte_m, gte_s] = display_glauber_data(fname, R)

numr = length(R);

[params, ~] = read_glauber_file(fname{1,1})
mi = zeros(numr, params.T);
te = zeros(numr, params.T);
gte = zeros(numr, params.T);
binarygte = zeros(numr, params.T);
reducedgte = zeros(numr, params.T);
reducedgtemanual = zeros(numr, params.T);
energy = zeros(numr, params.T);
magnetisation = zeros(numr, params.T);
interfacial_lengths = zeros(numr, params.T);
energy_hist = zeros(numr, params.L^2*2, params.T);

for r=1:numr
    [data] = calc_glauber_data(fname{1,r}, []);
    mi(r,:) = data.mi; te(r,:) = data.te; gte(r,:) = data.gte;
    binarygte(r,:) = data.binarygte; reducedgte(r,:) = data.reducedgte;
    magnetisation(r,:) = data.magnetisation;
    interfacial_lengths(r,:) = data.interfacial_lengths;
    energy_hist(r,:,:) = data.energy_hist;
end

Tc = 1/log(1+sqrt(params.numStates));
if params.pottsVersion == 'b'
	Tc = Tc * 2;
end
% For q=10
TcApprox = [0.70171, 0.70131, 0.70124];

mi_m = mean(mi,1)';
mi_s = std(mi,0,1)'/sqrt(numr);
te_m = mean(te,1)';
te_s = std(te,0,1)'/sqrt(numr);
gte_m = mean(gte,1)';
gte_s = std(gte,0,1)'/sqrt(numr);
binarygte_m = mean(binarygte,1)';
binarygte_s = std(binarygte,0,1)'/sqrt(numr);
reducedgte_m = mean(reducedgte,1)';
reducedgte_s = std(reducedgte,0,1)'/sqrt(numr);
reducedgtemanual_m = mean(reducedgtemanual,1)';
reducedgtemanual_s = std(reducedgtemanual,0,1)'/sqrt(numr);
energy_m = mean(energy,1)';
energy_s = std(energy,0,1)'/sqrt(numr);
magnetisation_m = mean(magnetisation,1)';
magnetisation_s = std(magnetisation,0,1)'/sqrt(numr);
ifl_m = mean(interfacial_lengths, 1)';
ifl_s = std(interfacial_lengths,0,1)'/sqrt(numr);


figure;
hold on;
ax1 = gca;
hold(ax1, 'on');
h=[];
grid('on');
tidx = params.T_vals>0;%>=Tc;
tvals = params.T_vals;

h = plot_error_bars(h, params.T_vals(tidx), mi_m(tidx), mi_s(tidx), '-ro', 'MI');
h = plot_error_bars(h, params.T_vals(tidx), te_m(tidx), te_s(tidx), '-go', 'TE');
h = plot_error_bars(h, params.T_vals(tidx), gte_m(tidx), gte_s(tidx), '-bo', 'GTE');
h = plot_error_bars(h, params.T_vals(tidx), binarygte_m(tidx), binarygte_s(tidx), '--rs', 'Binary GTE');
h = plot_error_bars(h, params.T_vals(tidx), reducedgte_m(tidx), reducedgte_s(tidx), '--bd', 'Reduced GTE');
h = plot_error_bars(h, params.T_vals(tidx), ifl_m(tidx), ifl_s(tidx), '--kx', 'Interfacial Length (ensemble average');
titleStr = sprintf('Simulation Metrics: q=%d, L=%d, U=%d, S=%d, T_c=%f, samples=%d', params.numStates, params.L, params.U, params.S, Tc, params.samples);
title(titleStr);
xlabel('T');
ylabel('Bits');

xlims =[min(params.T_vals),max(params.T_vals)];
line([Tc, Tc], ylim);
idx=1;
if params.L==64
    idx=1;
elseif params.L==128
    idx=2;
elseif params.L==512
    idx=3;
end
if params.numStates ==10
    line([TcApprox(idx), TcApprox(idx)], ylim,'LineStyle', '--');
end
xlim(xlims);

ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
hold(ax2, 'on');
xlim(ax2, xlims);
h(end+1) = errorbar(ax2, params.T_vals(tidx), magnetisation_m(tidx), magnetisation_s(tidx), '-ko', 'DisplayName', 'Magnetisation');
ylabel('Magnetisation');

legend(h, 'Location', 'Best');
hold off;

figure;
hold on;
E = [1:params.L^2*2];
EN = E *-1/params.L^2;
title('Energy histogram')
for t = 1:params.T
    plotStr = sprintf('T=%f', params.T_vals(t));
    plot(EN, energy_hist(1,:,t), 'DisplayName', plotStr);
end
hold off;

function h1 = plot_error_bars(h, x, y, ybar, linespec, name)
    h1 = h;
    if any(y)
        h1(end+1) = errorbar(x, y, ybar, linespec, 'DisplayName', name);
    end