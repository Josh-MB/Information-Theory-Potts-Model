function [data] = calc_glauber_data(fname, additional_runs)

[params, file_data] = read_glauber_file(fname);
total_runs=1;
% Used to create an ensemble
if any(additional_runs)
    for num1 = additional_runs
        total_runs = total_runs+1;
        slashes = strfind(fname, '/');
        tmpname = sprintf('%s%d%s', fname(1:slashes(end)-2), num1, fname(slashes(end):end));
        fprintf('Trying to open file %s\n', tmpname);
        [params_tmp, data_temp] = read_sim_file(tmpname);
        % Need to check params are equal
        file_data.mihist = file_data.mihist + data_temp.mihist;
        file_data.tehist = file_data.tehist + data_temp.tehist;
        file_data.gtehist = file_data.gtehist + data_temp.gtehist;
        file_data.binarygtehist = file_data.binarygtehist + data_temp.binarygtehist;
        file_data.reducedgtehist = file_data.reducedgtehist + data_temp.reducedgtehist;
        file_data.magnetisation = file_data.magnetisation + data_temp.magnetisation;
        file_data.energy_hist = file_data.energy_hist + data_temp.energy_hist;
    end
end
fprintf('\nHist made up of %d runs\n', total_runs);
data.magnetisation = file_data.magnetisation./total_runs;

data.tau_mean = file_data.tau_mean;
data.tau_std_error = file_data.tau_std_error;
data.energy_hist = file_data.energy_hist;
data.interfacial_lengths = file_data.interfacial_lengths;

data.mi = zeros(1,params.T);
data.te = zeros(1,params.T);
data.gte = zeros(1,params.T);
data.binarygte = zeros(1,params.T);
data.reducedgte = zeros(1,params.T);

for idx=1:params.T
    tmp = file_data.mihist(:,:,idx);
    P = tmp/sum(tmp(:));
    data.mi(idx) = entropy(sum(P,1)) + entropy(sum(P,2)) - entropy(P);
end

for idx=1:params.T
    % tmp is X',X,Y
    tmp = file_data.tehist(:,:,:,idx);
    P = tmp/sum(tmp(:));
    % PP is X',X
    PP = sum(P,3);
    data.te(idx) = centropy(PP,1) - centropy(P,1);
end

for idx=1:params.T
    tmp = file_data.gtehist(:,:,:,:,:,:,idx);
    P = tmp/sum(tmp(:));
    PP = sum(P(:,:,:),3);
    data.gte(idx) = centropy(PP,1) - centropy(P,1);
end

for idx=1:params.T
    tmp = file_data.binarygtehist(:,:,:,:,:,:,idx);
    P = tmp/sum(tmp(:));
    PP = sum(P(:,:,:),3);
    data.binarygte(idx) = centropy(PP,1) - centropy(P,1);
end

for idx=1:params.T
    % tmp is X',X,Delta
    tmp = file_data.reducedgtehist(:,:,:,idx);
    P = tmp/sum(tmp(:));
    % PP is X',X
    PP = sum(P,3);
    data.reducedgte(idx) = centropy(PP,1) - centropy(P,1);
end

data.mi = data.mi.*log2(exp(1));
data.te = data.te.*log2(exp(1));
data.gte = data.gte.*log2(exp(1));
data.reducedgte = data.reducedgte.*log2(exp(1));
data.binarygte = data.binarygte.*log2(exp(1));
