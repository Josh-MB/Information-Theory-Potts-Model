function [paramsGTE, tglT, T, gte, en, nlogp, enOrig, nlogpOrig, tglTgauss,tglTb, gte_binary, mag, mag_count, ple, magnetisation, ipw, mi, tglTs] = calc_dos_data(gteProfile, dosFile, draw, Eo, Ed)

fprintf('Reading: %s\n', gteProfile);
% Get DoS Data
[params, g, e, en] = readDoSFile(dosFile);

% Get GTE Profile Data
fid = fopen(gteProfile);

assert(check_magic(fid));
[paramsGTE.fileVer, c] = fread(fid,1,'int32'); assert(c==1);
[paramsGTE.L, c] = fread(fid,1,'uint64'); assert(c==1);
[paramsGTE.U, c] = fread(fid,1,'uint64'); assert(c==1);
[paramsGTE.T, c] = fread(fid,1,'double'); assert(c==1);
[paramsGTE.seed, c] = fread(fid, 1, 'uint64'); assert(c==1);
[paramsGTE.numStates, c] = fread(fid, 1, 'char'); assert(c==1);
[paramsGTE.samples, c] = fread(fid, 1, 'int32'); assert(c==1);
T=paramsGTE.T;
% GTE data
N=paramsGTE.L.^2;
numE1 = 32*32*2+1; % For scaling
%numE1 = 2*(paramsGTE.L^2)+1; % For not scaling
gte_binary = zeros(numE1, 1);
gte_binary_count = zeros(numE1,1);
[gte, c] = fread(fid, numE1, 'double'); assert(c==numE1);
[gteCounts, c] = fread(fid, numE1, 'uint64'); assert(c==numE1);
mi = zeros(numE1, 1);
miCounts = zeros(numE1, 1);
magnetisation = zeros(numE1, 1);
mag_count = zeros(numE1, 1);
gte_binary = zeros(numE1, 1);
gte_binary_count = zeros(numE1, 1);
gte_sweep = zeros(numE1, 1);
gte_sweep_count = zeros(numE1, 1);
    
if paramsGTE.fileVer >= 2
    [mi, c] = fread(fid, numE1, 'double'); assert(c==numE1);
    [miCounts, c] = fread(fid, numE1, 'uint64'); assert(c==numE1);
    [magnetisation, c] = fread(fid, numE1, 'double'); assert(c==numE1);
    [mag_count, c] = fread(fid, numE1, 'uint64'); assert(c==numE1);
    if paramsGTE.fileVer >= 3
        [gte_binary, c] = fread(fid, numE1, 'double'); assert(c==numE1);
        [gte_binary_count, c] = fread(fid, numE1, 'uint64'); assert(c==numE1);
        if paramsGTE.fileVer >= 4
            [gte_sweep, c] = fread(fid, numE1, 'double'); assert(c==numE1);
            [gte_sweep_count, c] = fread(fid, numE1, 'uint64'); assert(c==numE1);
        end
    end
end

fclose(fid);

logpOrig = g - e./paramsGTE.T;
enOrig = e ./ paramsGTE.L^2;
nlogpOrig = exp(logpOrig - max(logpOrig));

% If using a larger L value, have to reinsert the "missing values" near the
% end, then we average over the requisite number of bins, to collapse the 
% higher L down to L=32. Last value (ground state) does not get averaged
% and instead is appended to the end, then remove the missing values again
if params.L > 32
    numE2= 2*(params.L^2)+1;
    g2=zeros(numE2,1);
    g2(1:length(g)) = g;
    e2=zeros(numE2,1);
    e2(1:length(e)) = e;
    g2(numE2) = g2(numE2);
    g2(numE2-1) = g2(numE2);
    g2(numE2-2) = g2(numE2);
    g2(numE2-3) = g2(numE2);
    g2(numE2-5) = g2(numE2-4);
    e2(end) = e(end);
    e2(numE2-1) = e2(numE2)+1;
    e2(numE2-2) = e2(numE2)+2;
    e2(numE2-3) = e2(numE2)+3;
    e2(numE2-4) = e2(numE2)+4;
    e2(numE2-5) = e2(numE2)+5;
    nscale = (params.L / 32)^2;
    rescaledg = arrayfun(@(i) mean(g2(i:i+nscale-1)),1:nscale:length(g2)-nscale+1)';
    rescalede = arrayfun(@(i) mean(e2(i:i+nscale-1)),1:nscale:length(e2)-nscale+1)';
    rescaledg = [rescaledg; g2(end)];
    rescalede = [rescalede; e2(end)];
    g = rescaledg;
    e = rescalede;
    g(numE1-1,:) = [];
    g(numE1-2,:) = [];
    g(numE1-3,:) = [];
    g(numE1-5,:) = [];
    e(numE1-1,:) = [];
    e(numE1-2,:) = [];
    e(numE1-3,:) = [];
    e(numE1-5,:) = [];
    %e=floor(e);
    e1 = e./ (params.L^2 / 32^2);
    %e1=floor(e1);
    en = e./ (params.L^2);
    %en = e./(32*32);
    %en1 = e1./(32*32);
end

%params.L
%e(end-10:end)
%e1(end-10:end)
   
logp = g - e./paramsGTE.T;
nlogp = exp(logp - max(logp));

ple = double_gaussian(T, en, paramsGTE.numStates, params.L, 1, Eo, Ed);
ple=ple./max(ple);
%if params.L <= 32
    gte(numE1-1,:) = [];
    gte(numE1-2,:) = [];
    gte(numE1-3,:) = [];
    gte(numE1-5,:) = [];
    gteCounts(numE1-1) = [];
    gteCounts(numE1-2) = [];
    gteCounts(numE1-3) = [];
    gteCounts(numE1-5) = [];
    mi(numE1-1,:) = [];
    mi(numE1-2,:) = [];
    mi(numE1-3,:) = [];
    mi(numE1-5,:) = [];
    miCounts(numE1-1) = [];
    miCounts(numE1-2) = [];
    miCounts(numE1-3) = [];
    miCounts(numE1-5) = [];
    magnetisation(numE1-1,:) = [];
    magnetisation(numE1-2,:) = [];
    magnetisation(numE1-3,:) = [];
    magnetisation(numE1-5,:) = [];
    mag_count(numE1-1) = [];
    mag_count(numE1-2) = [];
    mag_count(numE1-3) = [];
    mag_count(numE1-5) = [];
    gte_binary(numE1-1,:) = [];
    gte_binary(numE1-2,:) = [];
    gte_binary(numE1-3,:) = [];
    gte_binary(numE1-5,:) = [];
    gte_binary_count(numE1-1) = [];
    gte_binary_count(numE1-2) = [];
    gte_binary_count(numE1-3) = [];
    gte_binary_count(numE1-5) = [];
    gte_sweep(numE1-1,:) = [];
    gte_sweep(numE1-2,:) = [];
    gte_sweep(numE1-3,:) = [];
    gte_sweep(numE1-5,:) = [];
    gte_sweep_count(numE1-1) = [];
    gte_sweep_count(numE1-2) = [];
    gte_sweep_count(numE1-3) = [];
    gte_sweep_count(numE1-5) = [];
    
%end

%gte(nlogp < 0.0001) = 0;
% Multiply gte(E,T) by exp term, and sum over E. Gives numerator
%gte=flipud(gte);
a = sum(gte .* nlogp);
% Sum exp term over E. Gives denominator
b = sum(nlogp);
tglT = a ./ b;
tglT = tglT * log2(exp(1));

tglTgauss = sum(gte.*ple)./sum(ple);
tglTgauss = tglTgauss * log2(exp(1));

ipw = sum(mi.*nlogp) ./ b;
ipw = ipw * log2(exp(1));

tglTb = sum(gte_binary.*nlogp) ./ b;
tglTb = tglTb * log2(exp(1));

tglTs = sum(gte_sweep.*nlogp) ./ b;
tglTs = tglTs * log2(exp(1));

magnetisation(isnan(magnetisation)) = 0;
magnetisation =  abs(magnetisation);

mag = sum(magnetisation .* nlogp) ./ b;

if draw
    fprintf('T=%f\n', T);
    figure;
    hold on;
    plot(en, gte, '-r');
    xlabel('E/N');
    ylabel('\sum_E Tgl(E,T)g(E)e^{-E/k_bT_c} / \sum_E g(E)e^{-E/k_bT_c}');
    title('GTE');
    hold off;
    
    figure;
    hold on;
    plot(en, gteCounts);
    xlabel('E/N');
    ylabel('Count');
    title('Energy histogram counts');
    hold off;
    
    figure;
    hold on;
    plot(en, nlogp, '-b', 'DisplayName', 'DoS');
    plot(en, ple, '-r', 'DisplayName', 'Double Gaussian');
    xlabel('E/N');
    ylabel('P(E,T)');
    title('Canonical Distribution');
    hold off;
end
