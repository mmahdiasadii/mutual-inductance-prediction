%% Random Numbers
tic
clear; clc;
DataNum=15000;
n1 = randi([4, 10], DataNum, 1);
n2 = randi([4, 10], DataNum, 1);
nt1 = randi([5, 15], DataNum, 1);
nt2 = randi([5, 15], DataNum, 1);
s1 = Random_Numbers(1e-3, 3e-3, DataNum, 1);
s2 = Random_Numbers(1e-3, 3e-3, DataNum, 1);

ri1 = Random_Numbers(30e-3, 80e-3, DataNum, 1);
ri2 = Random_Numbers(30e-3, 80e-3, DataNum, 1);
% 
ro1 = ri1+(nt1-1).*s1./cos(pi./n1);
% ro2 = ro1.*Random_Numbers(0.8, 1, 1, 1);

% s2 = (ro2 - ri2).*cos(pi./n2)./(nt2-1);

phi = Random_Numbers(-pi/12, pi/12, DataNum, 1);
tetha = Random_Numbers(-pi/12, pi/12, DataNum, 1);
psi = Random_Numbers(-2*pi./n2, 2*pi./n2, DataNum, 1);
l = Random_Numbers(-ro1*sqrt(2)/2, ro1*sqrt(2)/2, DataNum, 1);
w = Random_Numbers(-ro1*sqrt(2)/2, ro1*sqrt(2)/2, DataNum, 1);
h = Random_Numbers(50e-3, 120e-3, DataNum, 1);

out_mat = [n1, n2, nt1, nt2, s1, s2, ri1, ri2,...
    phi, tetha, psi, l, w, h];

%% Mutual Inductance Calculations
M = zeros(DataNum, 1);
for i=1:DataNum
    M(i, 1) = Mutual_Inductance(n1(i), n2(i), nt1(i),...
        nt2(i), s1(i), s2(i), ri1(i), ri2(i), phi(i), tetha(i),...
        psi(i), l(i), w(i), h(i));
end
%% CSV file write

columnLabels = {'n1', 'n2', 'nt1', 'nt2', 's1', 's2', 'ri1', 'ri2',...
    'phi', 'tetha', 'psi', 'l', 'w', 'h', 'M'};
dataTable = table(out_mat(:,1), out_mat(:,2), out_mat(:,3),...
    out_mat(:,4), out_mat(:,5), out_mat(:,6), out_mat(:,7),...
    out_mat(:,8), out_mat(:,9), out_mat(:,10), out_mat(:,11),...
    out_mat(:,12), out_mat(:,13), out_mat(:,14), M(:,1),...
    'VariableNames', columnLabels);
fileName = 'output.csv';
writetable(dataTable, fileName);
disp(['Data has been written to ' fileName]);
toc