function [conclusion] = Catalogue_2_Summary( file )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

daten = importdata(file);
eq = unique(daten.data(:,1));

conclusion = [];
for i = 1:1:length(eq)
    daten_red = daten.data(daten.data(:,1) == eq(i),: );
    IS = daten_red(:,10);
    IX = daten_red(:,2);
    ID = daten_red(:,1);
    


    for k = 1:1:12
    summ(k) = sum(floor(IS) == k);
    end
    
    total = sum(summ);
 
 conclusion = [conclusion; eq(i) IX(1) summ total ];
end
filename =  strcat({ 'summary-'},{file});
 fileID = fopen(filename{1},'w');
 fprintf(fileID,'%6s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n',...
     'EQ ID', 'Ix','1','2','3','4','5','6','7','8','9','10','11','12', 'Total');
 fprintf(fileID,'%6d %12.4f %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d\n', conclusion');
 fclose(fileID);
end

