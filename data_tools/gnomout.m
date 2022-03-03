function [data_original, datanew, Pairdist] = gnomout(file)
%function [data_original, datanew, Pairdist] = gnomout(file)

if ~exist('file')
	[file, datadir]=uigetfile('*.out','Select data file');
	if file==0 return; end
    file = [datadir, file];
end

fid = fopen(file,'r');

firstline = ' ';
startn = ' ';
data  = [];
data2 = [];
numofrun = 1;

while startn ~= 'OK'
   firstline = fgetl(fid);
   if (length(firstline) >= 7)&(firstline(1:7) == '      S')
       firstline = fgetl(fid);
       firstline = fgetl(fid);
       break
   end
end       

i = 1;
k = 1;

while isempty(firstline) ~= 1
    %firstline = retrival_num(firstline);
    dataline = sscanf(firstline, '%f');
    if length(dataline) == 5
        data(k, 1) = dataline(1);
        data(k, 2) = dataline(2);
        data(k, 3) = dataline(3);
        data(k, 4) = dataline(4);
        data(k, 5) = dataline(5);
        k = k+1;
    elseif length(dataline) == 2
        data2(i, 1) = dataline(1);
        data2(i, 2) = dataline(2);
        i = i+1;
    end
        firstline = fgetl(fid);
end

data_original = [data(:, 1), data(:,2)];
datanew(:,1) = [data2(:,1)
    data(:,1)];
datanew(:,2) = [data2(:,2)
    data(:,5)];

clear data
while startn ~= 'OK'
   firstline = fgetl(fid);
   if (length(firstline) >= 8)&(firstline(1:8) == '       R')
       firstline = fgetl(fid);
       firstline = fgetl(fid);
       break
   end
end       

k = 1;

while (isempty(firstline) ~= 1) & (isstr(firstline) ~= 0)
    %firstline = retrival_num(firstline);
    dataline = sscanf(firstline, '%f');
    
    if length(dataline) == 3
        data(k, 1) = dataline(1);
        data(k, 2) = dataline(2);
        data(k, 3) = dataline(3);
        k = k+1;
    end
    firstline = fgetl(fid);
end

Pairdist = data;

fclose(fid);