function good_magic = check_magic(fid,verbose)

good_magic = false;

if nargin < 2 || isempty(verbose), verbose = false; end

tline = '';
bmag = false;
while ischar(tline)
    tline = fgetl(fid);
    if verbose, disp(tline); end
    if strcmp(tline,'BEGIN MAGIC')
        bmag = true;
        break;
    end
end
if ~bmag, error('magic does not begin'); end

[xc, c] = fread(fid,3,'char');   if ~(c == 3 && isequal(xc,['b';'a';'h'])), fprintf(2,'WARNING - bad magic: char\n');    return; end
[xs, c] = fread(fid,1,'int64');  if ~(c == 1 && xs == 12345678987654321),   fprintf(2,'WARNING - bad magic: size_t\n');  return; end
%[x1, c] = fread(fid,1,'int16');  if ~(c == 1 && x1 == 1234),                fprintf(2,'WARNING - bad magic: short\n');   return; end
%[y1, c] = fread(fid,1,'int16');  if ~(c == 1 && y1 == -x1),                 fprintf(2,'WARNING - bad magic: -short\n');  return; end
[x2, c] = fread(fid,1,'int32');  if ~(c == 1 && x2 == 12345678),            fprintf(2,'WARNING - bad magic: int\n');     return; end
[y2, c] = fread(fid,1,'int32');  if ~(c == 1 && y2 == -x2),                 fprintf(2,'WARNING - bad magic: -int\n');    return; end
[x3, c] = fread(fid,1,'int64');  if ~(c == 1 && x3 == 1234567898765432),    fprintf(2,'WARNING - bad magic: long\n');    return; end
[y3, c] = fread(fid,1,'int64');  if ~(c == 1 && y3 == -x3),                 fprintf(2,'WARNING - bad magic: -long\n');   return; end
[x4, c] = fread(fid,1,'double'); if ~(c == 1 && x4 == 12345678.987654321),  fprintf(2,'WARNING - bad magic: double\n');  return; end
[y4, c] = fread(fid,1,'double'); if ~(c == 1 && y4 == -x4),                 fprintf(2,'WARNING - bad magic: -double\n'); return; end

if verbose, fprintf('\t<magic ok>\n'); end

tline = fgetl(fid);
if verbose, disp(tline); end
assert(ischar(tline),'magic gap');

tline = fgetl(fid);
if verbose, disp(tline); end
assert(ischar(tline),'magic does not end');
assert(strcmp(tline,'END MAGIC'),'magic ends badly');

good_magic = true;
