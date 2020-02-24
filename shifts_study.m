% This is 1/3 part of analysing toolkit, used for ploting numbers of drifts inside each trail;
% Usage: input `number_1`[array]
%              `trails_1`[array], but numbers should be equivalent to trails
% Author: Qi Wang, AGYU, Max Planck Institute für Biological Kybernetics,
% 72076 Tübingen,Germany


% baseline drifts

% numbers_1 = [0,0,0,0,3,0,0,3,0,4,0,2,3,0,3,0,0,2,1,2,1];%0523
numbers_1 = [0,0,0,0,0,0,0,0,1,1,0,0,3,1,2,2,1,1,0,2,3,1,0,0,2,1,2,2,0,0,0];%0530
trails_1 = linspace(1,length(numbers_1),length(numbers_1));
% num of shifts each trail
p = plot(trails_1,numbers_1,'--','LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[0.8500 0.3250 0.0980])
p.Marker = 'o';
xlim([0,21]);
ylim([0,5,]);
yticks([0,1,2,3,4]);
xlabel('Trails');
ylabel('Quantity of shifts');
grid on;
legend('Baseline shifts');

% mean xc coeff of each trail
m_1 = mean(pks(:,:,3),1)% only see 3rd slice
p_1 = plot(trails_1,m_1,'--','LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[0.8500 0.3250 0.0980])
p_1.Marker = 'o';
xlim([0,21]);
yticks([0,1,2,3,4]);
xlabel('Trails');
ylabel('XC coeff mean');
grid on;
legend('mean XC coeff');