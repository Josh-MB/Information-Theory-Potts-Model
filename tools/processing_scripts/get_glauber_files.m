function [fname,R] = get_glauber_files(subdir)
% Grab all .bin files in specified subdirectory and sort
% according to runs

fdir = fullfile(getenv('DATADIR'),'.',subdir);

files = fullfile(fdir,'*.bin');
fprintf('finding files : %s\n\n',files);
fname = struct2cell(dir(files));
fname = fname(1,:)';
nfiles = length(fname)

R = zeros(nfiles,1);
for f = 1:nfiles
    sfname = fname{f};
    i(1:2) = strfind(sfname,'_');
    i(3) = strfind(sfname,'.bin');
    %E(f) = str2double(sfname(i(1)+1:i(2)-1));
    R(f) = str2double(sfname(i(2)+1:i(3)-1));
    fname{f} = fullfile(fdir,fname{f});
end
assert(min(R) == 1);
numr = max(R);
nume = nfiles/numr;
assert(nume*numr == nfiles);
fname = reshape(fname,numr,nume)';

R = reshape(R,numr,nume);

R = R(:,1);
