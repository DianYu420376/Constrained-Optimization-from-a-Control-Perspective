data1 = load('case39_history.mat');
data2 = load('case39_history_momentum.mat');
data3 = load('case39_history_newton.mat');
line1 = data1.KKT_gaps;
line2 = data2.KKT_gaps;
line3 = data3.KKT_gaps;
plot(log10([line1, line2,line3]))