function [ catalogue ] = DBMI04_2_Catalogue(files)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



catalogue = [];
for i = 1:1:length(files)
a = importdata(files{i});

epicentro = cell2mat(a(4,2:3));
eq_info = cell2mat(a(6,1));

epi_lat = str2double(epicentro(12:17));
epi_long = str2double(epicentro(20:25));
mw = str2double(epicentro(35:38));
mw_err = str2double(epicentro(41:44));
nobs = str2double(eq_info(4:5));


% convert Ix to numeric. if uncertain assign mean if string set NaN
ix = eq_info(14:end);
ix = ix(~isspace(ix));
ix = ix(~isletter(ix));
Ix = str2double(ix);
Ix(strcmp(ix, '1-2'))=1.5;
 Ix(strcmp(ix, '2-3'))=2.5;
 Ix(strcmp(ix, '3-4'))=3.5;
 Ix(strcmp(ix, '4-5'))=4.5;
 Ix(strcmp(ix, '5-6'))=5.5;
 Ix(strcmp(ix, '6-7'))=6.5;
 Ix(strcmp(ix, '7-8'))=7.5;
 Ix(strcmp(ix, '8-9'))=8.5;
 Ix(strcmp(ix, '9-10'))=9.5;
 Ix(strcmp(ix, '10-11'))=10.5;
 Ix(strcmp(ix, '11-12'))=11.5;


data = a(9:end,3:5);
lat = cell2mat(cellfun(@str2num, data(:,1),'UniformOutput', false));
lon = cell2mat(cellfun(@str2num, data(:,2),'UniformOutput', false));
is = data(:,3);

% convert Is to numeric. if uncertain assign mean. if string set NaN
Is = str2double(is);
 Is(strcmp(is, '1-2'))=1.5;
 Is(strcmp(is, '2-3'))=2.5;
 Is(strcmp(is, '3-4'))=3.5;
 Is(strcmp(is, '4-5'))=4.5;
 Is(strcmp(is, '5-6'))=5.5;
 Is(strcmp(is, '6-7'))=6.5;
 Is(strcmp(is, '7-8'))=7.5;
 Is(strcmp(is, '8-9'))=8.5;
 Is(strcmp(is, '9-10'))=9.5;
 Is(strcmp(is, '10-11'))=10.5;
 Is(strcmp(is, '11-12'))=11.5;
 
 IS = Is(~strcmp(data(:,1),''));
 
 
result =[lat, lon, IS, repmat([ Ix,epi_lat, epi_long, mw, mw_err,i], ...
    length(lat),1)];
result = result(~isnan(result(:,3)),:);

clear lat lon IS Ix epi_lat epi_long mw mw_err i

catalogue = [catalogue; result];
end

end

