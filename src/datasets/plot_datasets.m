%% Import datasets as Matlab table

confirmed = readtable('confirmed_cases_metr_france.csv', 'ReadVariableNames', true);
death = readtable('death_cases_metr_france.csv', 'ReadVariableNames', true);
recovered = readtable('recovered_cases_metr_france.csv', 'ReadVariableNames', true);

%% Plot datasets on a single axes

figure(1); clf;
semilogy(confirmed.Date, confirmed.Value,'LineWidth',2)
hold on % do not replace data on plot
semilogy(death.Date, death.Value,'LineWidth',2)
semilogy(recovered.Date, recovered.Value,'LineWidth',2)

% text(datetime('2020-01-29'),1e3,'https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases')

% add gov measures
line(datetime({'2020-03-16','2020-03-16'}),[1e0,1e5],'Color','r', 'LineStyle', '-.');
line(datetime({'2020-03-17','2020-03-17'}),[1e0,1e5],'Color','b', 'LineStyle', '--');

legend('Confirmed', 'Deaths', 'Recovered','Schools closed','Isolation','Location','northwest')
xlabel('Date')
ylabel('Number')
title({'Covid-19 in France';'source: https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases'})
