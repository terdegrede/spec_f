clc, clearvars, close all

load pop.dat
load vmat.dat

[M, N] = size(pop);
P = pop;
B = vmat;


dB = 0.05*B;     % Treshold for reproductive isolation
dS = dB*N/B;

DR = pdist2(P, P, 'hamming')*size(P,2); % Full genome

figure('WindowState','maximized', 'Visible','on');
axes('Units', 'normalized', 'Position', [0 0 1 1]);
%% dS
% To form clusters (Blue)
DR(DR>dS) = 0;
Gr2 = graph(DR,'upper');
axes('Position',[-0.25 0 1 1]);
GR = plot(Gr2);
axis square
GR.NodeLabel = {};
GR.EdgeAlpha = 1;
GR.EdgeColor = [0 0 1]; % blue
xdata = GR.XData;
ydata = GR.YData;

%% dB
% To link reproductive compatible (Red)
PB = P(:, 1:B);

DRB = pdist2(PB, PB, 'hamming')*B; % Mating trait
DRB(DRB>dB) = 0;
GrB = graph(DRB,'upper');
axes('Position',[0.25 0 1 1]);
GRB = plot(GrB);
axis square
GRB.NodeColor = 'black';
GRB.NodeLabel = {};
GRB.EdgeAlpha = 1;
GRB.XData = xdata;
GRB.YData = ydata;
GRB.LineStyle = '--';
GRB.EdgeColor = [1 0 0]; % red