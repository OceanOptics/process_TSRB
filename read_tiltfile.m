function [data,datatest,tilt_tt2] = read_tiltfile(filename,filedir,tiltlim)

filename=[filedir,filename];

% Read file header lines to get number of data columns
fid = fopen(filename);
A = textscan(fid,'%s',5,'delimiter','\n');
Acell = A{1};
colnames = Acell(4);
findcommas = regexp(colnames,',');
numcols = length(findcommas{1})-1;

number_of_string_columns = 2;

formatpattern = [];

for ii = 1:numcols
    formatpattern = [formatpattern '%f'];
end
for ii = 1:number_of_string_columns
    formatpattern = [formatpattern '%s'];
end

% read in the file contents
A = textscan(fid,formatpattern,'headerlines',5,'delimiter',',');
fclose(fid);

%find the idices where the tilt in both x and y directions is less than
%user defined number of degrees
xcol = regexp(colnames,',TILT(X');
xcol_ind = (find(xcol{1} == findcommas{1}))+1;
ycol = regexp(colnames,',TILT(Y');
ycol_ind = (find(ycol{1} == findcommas{1}))+1;

tiltx=A{xcol_ind};
tilty=A{ycol_ind};

keeptilt=[];
for itilt=1:(max(size(tiltx)))
    if (tiltx(itilt)<tiltlim) && (tiltx(itilt)>-1*tiltlim)
        if (tilty(itilt)<tiltlim) && (tilty(itilt)>-1*tiltlim)
            keeptilt=[keeptilt,itilt];
        end
    end
end

data=tiltx(keeptilt);
datatest=tilty(keeptilt);

% read date and time data (renamed as 'tt2') as formatted by SatCon
for i = 1:number_of_string_columns
    tt2_temp(:,i) = A{1,numcols+i};
end

% % Convert SatCon date/tt2 to Matlab datenumber format
for i = 1:size(tt2_temp,1)
    B = [sscanf(tt2_temp{i,1},'%f-%f')' sscanf(tt2_temp{i,2},'%f:%f:%f')'];
    C = datevec(datenum([B(1) 0 B(2)]));
    C(4:6) = [B(3) B(4) B(5)];
    tilt_tt2_all(i,1) = datenum(C);
end

%get only the times with good tilt values
tilt_tt2=tilt_tt2_all(keeptilt);

end
