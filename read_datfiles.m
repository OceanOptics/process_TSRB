function [data,wl,tt2,SN,INTtt2,filetype] = read_datfiles(fnamdark,fnamlight,filedir,shutter)
% ------------------------------------------------------------------------
% Read data into Matlab.  Data files are assumed to have already been
% converted to comma-delimited text and calibration and immersion
% coefficients applied in SatCon.  Five header lines are assumed (to ensure
% this, CHECK "extra output header information" in SatCon parameters),
% and tt2 tags are assumed to have been processed and formatted (should be
% 2 last columns of .dat file).
%
% INPUTS:
% fnamdark  = string containing the path/filename of the dark-frame
%             data file for a single sensor
% fnamlight = string containing path/filename of the corresponding
%             light-frame data file
% filedir   = directory of .dat files
% shutter   = 1 if the Lu shutter was not working
%
% OUTPUTS:
% data      = radiometric data with tt2-interpolated dark-frame
%             values subtracted
% wl        = array of wavelengths associated with radiometric data
% tt2       = Matlab datenumber-formated tt2 vector (sensor-
%             specific)
% SN        = Sensor serial number
% INTtt2    = vector of measurement integration tt2s (seconds)
% filetype  = type of radiometric measurement (Ed or Lu)
% ------------------------------------------------------------------------
fnamdark=[filedir,fnamdark];
fnamlight=[filedir,fnamlight];

if nargin == 3;shutter = 0;end

% Determine the number of columns in the file
fid = fopen(fnamdark);
line = fgetl(fid);line = fgetl(fid);line = fgetl(fid);line = fgetl(fid);
fclose(fid);
tmp=textscan(line, '%s','delimiter',',');
allcols = length(tmp{:});

% 1) READ IN THE DARK-FRAME DATA FILE ...
fid = fopen(fnamdark);
formatpattern = [];
number_of_floating_point_columns = allcols-2;
number_of_string_columns= 2;            % inspect data file to determine
for ii = 1:number_of_floating_point_columns
    formatpattern = [formatpattern '%f'];
end
for ii = 1:number_of_string_columns
    formatpattern = [formatpattern '%s'];
end

% read in the file contents
A = textscan(fid,formatpattern,'headerlines',5,'delimiter',',');
fclose(fid);
% read numeric data (not all are output; see "rawheader" variables in
% OCRstuff.mat)
% if the file is empty, then set output as NaNs, otherwise write in the
% dark data and date/time info, then proceed to the light data file
if isempty(A{1})==1
    data=NaN;
    tt2=NaN;
    SN=NaN;
    INTtt2=NaN;
else
    
    % read in the fourth row of the data file to get header/wavelength information
    header=importdata(fnamdark,' ',4);
    head=header{3};
    starter=strfind(head,'3');
    xarr=head(starter(1):end);
    filetype=xarr(9:10);
    
    newx=strrep(xarr,filetype,'');
    newx2=strrep(newx,'),(',' ');
    ender=strfind(newx2,')');
    wvns=newx2(1:ender(1)-1);
    wldark=str2num(wvns);
    
    % write the dark file numerical data to a 'darkdata' array
    for i = 1:number_of_floating_point_columns
        darkdata(:,i) = A{1,i};
    end
    
    % read date and time data (renamed as 'tt2') as formatted by SatCon
    for i = 1:number_of_string_columns
        tt2(:,i) = A{1,number_of_floating_point_columns+i};
    end
    
    % Convert SatCon date/tt2 to Matlab datenumber format
    for i = 1:size(tt2,1)
        B = [sscanf(tt2{i,1},'%f-%f')' sscanf(tt2{i,2},'%f:%f:%f')'];
        C = datevec(datenum([B(1) 0 B(2)]));
        C(4:6) = [B(3) B(4) B(5)];
        darkdate(i,1) = datenum(C);
    end
    
    % clean up
    clear tt2 A B C;
    
    %repeat all of the above for light-frame datafile
    
    % fnamlight=HSE;
    fid = fopen(fnamlight);
    A = textscan(fid,formatpattern,'headerlines',5,'delimiter',',');
    fclose(fid);
end

if isempty(A{1})==1
    data=NaN;
    tt2=NaN;
    SN=NaN;
    INTtt2=NaN;
else
    
    % read in the fourth row of the data file to get header/wvn information
    header=importdata(fnamlight,' ',4);
    head=header{3};
    starter=strfind(head,'3');
    xarr=head(starter(1):end);
    filetype=xarr(9:10);
    
    newx=strrep(xarr,filetype,'');
    newx2=strrep(newx,'),(',' ');
    ender=strfind(newx2,')');
    wvns=newx2(1:ender(1)-1);
    wl=str2num(wvns);
    
    if wldark ~= wl
        error('wavelength arrays in dark and light file do not match')
    end
    
    for i = 1:number_of_floating_point_columns;
        lightdata(:,i) = A{1,i};
    end;
    
    for i = 1:2;
        tt2(:,i) = A{1,number_of_floating_point_columns+i};
    end;
    
    for i = 1:size(tt2,1);
        B = [sscanf(tt2{i,1},'%f-%f')' sscanf(tt2{i,2},'%f:%f:%f')'];
        C = datevec(datenum([B(1) 0 B(2)]));
        C(4:6) = [B(3) B(4) B(5)];
        lightdate(i,1) = datenum(C);
    end;
    % clean up
    clear tt2 A B C;
    
    % retrieve the serial number of the OCR sensor:
    SN = lightdata(1,2);
    
    % Prior to interpolating dark-frame data to light-frame tt2scale,
    % discard light frames collected before the first or after the last
    % dark frame. IF the dark data is coming from a separate data
    % collection (from when the Lu shutter was not working), force the dark
    % data to be on the time resolution of the light data.
    if shutter == 0;
        first = darkdate(1,1);
        last = darkdate(length(darkdate),1);
        good = intersect(find(lightdate>first),find(lightdate<last));
        
        % retrieve the vector of measurement integration tt2s (seconds)
        INTtt2 = lightdata(good,3);
        
        % retrieve the data. user must specify the columns where the
        % radiometer data are located:
        start_col = 5; % first colunm of radiometric data
        end_col = allcols-8; % last column of radiometric data
        lightdata = lightdata(:,start_col:end_col);
        darkdata = darkdata(:,start_col:end_col);
        
        % loop through each of the data columns and interpolate dark-frame
        % data to light-frame tt2scale
        for i = 1:size(lightdata,2);
            darkfill(:,i) = ...
                interp1(darkdate,darkdata(:,i),lightdate(good),'linear');
        end
        data = lightdata(good,:)-darkfill;
        tt2 = lightdate(good);
        
    else
        INTtt2 = lightdata(:,3);
        tt2 = lightdate;
        % retrieve the data. user must specify the columns where the
        % radiometer data are located:
        start_col = 5; % first colunm of radiometric data
        end_col = allcols-8; % last column of radiometric data
        lightdata = lightdata(:,start_col:end_col);
        darkdata = darkdata(:,start_col:end_col);
        
        % use the mean of all the dark data
        for iii = 1:length(darkdata(1,:))
            meandark(iii) = mean(darkdata(:,iii));
        end
        
        for ispec = 1:length(lightdata(:,1))
            data(ispec,:) = lightdata(ispec,:) - meandark;
        end
        
    end
    
    
end