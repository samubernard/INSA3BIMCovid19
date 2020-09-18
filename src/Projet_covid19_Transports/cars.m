opts = detectImportOptions('tmja-2017.csv');
M = readtable('tmja-2017.csv',opts);
plot(M.Var9,M.Var10,'.');
%Moyenne Autoroute A1
rows = (strcmp(M.(2),'A0001'));
t = table2array(M(rows,23));
moyA1 = mean(t);
%Moyenne autoroute A4
rows = (strcmp(M.(2),'A0004'));
t = table2array(M(rows,23));
moyA4 = mean(t);
%Moyenne autoroute A5
rows = (strcmp(M.(2),'A0005'));
t = table2array(M(rows,23));
moyA5 = mean(t);
%Moyenne autoroute A6
rows = (strcmp(M.(2),'A0006'));
t = table2array(M(rows,23));
moyA6 = mean(t);
%Moyenne autoroute A7
rows = (strcmp(M.(2),'A0007'));
t = table2array(M(rows,23));
moyA7 = mean(t);
%Moyenne autoroute A8
rows = (strcmp(M.(2),'A0008'));
t = table2array(M(rows,23));
moyA8 = mean(t);
%Moyenne autoroute A9
rows = (strcmp(M.(2),'A0009'));
t = table2array(M(rows,23));
moyA9 = mean(t);
%Moyenne autoroute A10
rows = (strcmp(M.(2),'A0010'));
t = table2array(M(rows,23));
moyA10 = mean(t);
%Moyenne autoroute A11
rows = (strcmp(M.(2),'A0011'));
t = table2array(M(rows,23));
moyA11 = mean(t);
%Moyenne autoroute A13
rows = (strcmp(M.(2),'A0013'));
t = table2array(M(rows,23));
moyA13 = mean(t);
%Moyenne autoroute A16
rows = (strcmp(M.(2),'A0016'));
t = table2array(M(rows,23));
moyA16 = mean(t);
%Moyenne autoroute A20
rows = (strcmp(M.(2),'A0020'));
t = table2array(M(rows,23));
moyA20 = mean(t);
%Moyenne autoroute A89
rows = (strcmp(M.(2),'A0089'));
t = table2array(M(rows,23));
moyA89 = mean(t);
%Sauvegarde dans un fichier csv
D = [1 4 5 6 7 8 9 10 11 13 16 20 89 ; 
    moyA1 moyA4 moyA5 moyA6 moyA7 moyA8 moyA9 moyA10 moyA11 moyA13 moyA16 moyA20 moyA89; 
    0.89*moyA1 0.89*moyA4 0.89*moyA5 0.89*moyA6 0.89*moyA7 0.89*moyA8 0.89*moyA9 0.89*moyA10 0.89*moyA11 0.89*moyA13 0.89*moyA16 0.89*moyA20 0.89*moyA89 ];
csvwrite('autoroutes.csv',D) ;