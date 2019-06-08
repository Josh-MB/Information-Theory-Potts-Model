function [params, data] = read_glauber_file(fname)

fprintf('Reading %s\n', fname);
fid = fopen(fname);

assert(check_magic(fid));
[params.L, c] = fread(fid,1,'uint64'); assert(c==1);
[params.U, c] = fread(fid,1,'uint64'); assert(c==1);
[params.S, c] = fread(fid,1,'uint64'); assert(c==1);
[params.samples,c] = fread(fid, 1, 'uint64'); assert(c==1);
[params.imode, c] = fread(fid,1,'int32'); assert(c==1);
[params.seed, c] = fread(fid,1,'uint64'); assert(c==1);
[params.numStates, c] = fread(fid,1,'char'); assert(c==1);
[params.pottsVersion, c] = fread(fid,1,'char'); assert(c==1);
[params.regime, c] = fread(fid,1,'int32'); assert(c==1);
%[params.orderedProp, c] = fread(fid, 1, 'double'); assert(c==1);

[params.T, c] = fread(fid,1,'int32'); assert(c==1);

[params.T_vals, c] = fread(fid, params.T, 'double'); assert(c==params.T);

[mihist,c] = fread(fid, params.numStates .^ 2 * params.T, 'uint64'); assert(c==params.numStates.^2*params.T);
data.mihist = reshape(mihist, [params.numStates, params.numStates, params.T]);

%[tehist,c] = fread(fid, params.numStates .^ 3 * params.T, 'uint64'); assert(c==params.numStates.^3*params.T);
%tehist = reshape(tehist, [params.numStates, params.numStates, params.numStates, params.T]);
data.tehist = zeros(params.numStates, params.numStates, params.numStates, params.T);

[gtehist,c] = fread(fid, params.numStates .^ 6 * params.T, 'uint64'); assert(c==params.numStates.^6*params.T);
gtehistReshape = ones(1,7) * params.numStates; gtehistReshape(7) = params.T;
data.gtehist = reshape(gtehist, gtehistReshape);

%[binarygtehist,c] = fread(fid, params.numStates .^ 2 * 16 * params.T, 'uint64'); assert(c==params.numStates.^2*16*params.T);
binarygtehistReshape = ones(1,7) * params.numStates; binarygtehistReshape(3:6)=2; binarygtehistReshape(7) = params.T;
%binarygtehist = reshape(binarygtehist, binarygtehistReshape);
data.binarygtehist = zeros(binarygtehistReshape);

[reducedgtehist,c] = fread(fid, params.numStates .^ 2 * 5 * params.T, 'uint64'); assert(c==params.numStates.^2*5*params.T);
%reducedgtehist = zeros(1,params.numStates.^2*9*params.T);
data.reducedgtehist = reshape(reducedgtehist, [params.numStates params.numStates 5 params.T]);

[data.magnetisation,c] = fread(fid, params.T, 'double'); assert(c==params.T);
[data.interfacial_lengths,c] = fread(fid, params.T, 'double'); assert(c==params.T);
%interfacial_lengths=zeros(1,params.T);
%magnetisation = zeros(1, params.T);

[energy_hist,c] = fread(fid, params.L^2*2*params.T, 'uint64'); assert(c==params.L^2*2*params.T);
data.energy_hist = reshape(energy_hist, [params.L^2*2, params.T]);

[data.tau_mean,c] = fread(fid, params.T, 'uint64'); assert(c==params.T);
[data.tau_std_error,c] = fread(fid, params.T, 'uint64'); assert(c==params.T);

fclose(fid);
