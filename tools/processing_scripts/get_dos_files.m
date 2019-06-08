%Get files from scaled output runs
function [fname,R,T] = get_dos_files(subdir)

fdir = fullfile(getenv('DATADIR'),'.',subdir);

files = fullfile(fdir,'*.bin');
fprintf('finding files : %s\n\n',files);
fname = struct2cell(dir(files));
fname = fname(1,:)';
nfiles = length(fname);

R = zeros(nfiles,1);
for f = 1:nfiles
    sfname = fname{f};
    i(1:4) = strfind(sfname,'_');
    i(5) = strfind(sfname,'.bin');
    %E(f) = str2double(sfname(i(1)+1:i(2)-1));
    R(f) = str2double(sfname(i(4)+1:i(5)-1));
    T(f) = str2double(sfname(i(3)+2:i(4)-1));
    fname{f} = fullfile(fdir,fname{f});
end
assert(min(R) == 1);
numr = max(R);
nume = nfiles/numr;
assert(nume*numr == nfiles);
fname = reshape(fname,numr,nume)';

R = reshape(R,numr,nume);
T = reshape(T,numr,nume);

R = R(:,1);
T = T(1,:);
