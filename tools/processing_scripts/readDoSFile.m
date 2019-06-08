function [ params, logg, e, en ] = readDoSFile( dosFile )
%READDOSFILE Reads in a DoS File
%   Reads in a DoS file and returns all data. Note there are multiple
%   versions of DoS file format, which are handled automatically by this
%   script.

fid = fopen(dosFile);

assert(check_magic(fid));
[dosFileVer, c] = fread(fid,1,'int32'); assert(c==1);
[params.L, c] = fread(fid,1,'uint64'); assert(c==1);
[params.U, c] = fread(fid,1,'uint64'); assert(c==1);
[params.seed, c] = fread(fid,1,'uint64'); assert(c==1);
if dosFileVer > 1
    [params.len, c] = fread(fid,1,'uint64'); assert(c==1);
    [params.rand_engine,c] = fread(fid,params.len,'char'); assert(c==params.len);
end
[params.numStates, c] = fread(fid,1,'char'); assert(c==1);
if dosFileVer < 3
	[params.Toffset, c] = fread(fid,1,'int32'); assert(c==1);
	[params.Tmin, c] = fread(fid, 1, 'double'); assert(c==1);
	[params.Tmax, c] = fread(fid, 1, 'double'); assert(c==1);
end
if dosFileVer > 1
    [params.factor,c] = fread(fid, 1, 'double'); assert(c==1);
    [params.iteration,c] = fread(fid,1,'uint64'); assert(c==1);
end
[params.gSize, c] = fread(fid, 1, 'uint64'); assert(c==1);

% Log of DoS
[logg, c] = fread(fid, params.gSize, 'double'); assert(c==params.gSize);
% Last histogram
[h, c] = fread(fid, params.gSize, 'uint64'); assert(c==params.gSize);
% Energy values
[e, c] = fread(fid, params.gSize, 'uint64'); assert(c==params.gSize);
e = -e/2;
en = e./(params.L*params.L);
fclose(fid);

end

