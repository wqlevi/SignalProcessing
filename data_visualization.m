% This is 2/3 part of analysing toolkit, used for ploting numbers of drifts inside each trail;
% Usage: 
% import pks and locs, which are animal specific XC coeff and lag time.
% Author: Qi Wang, AGYU, Max Planck Institute für Biological Kybernetics,
% 72076 Tübingen,Germany


% data import
load pks locs
nslice = 3; % 1 | 2 | 3 for 3 slices
for ii = 1:size(locs,2)
test_1((ii-1)*40+1:(ii-1)*40+40,1) = locs(:,ii,nslice);%lag times
test_1((ii-1)*40+1:(ii-1)*40+40,2) = pks(:,ii,nslice);%coeff
test_1((ii-1)*40+1:(ii-1)*40+40,3) = linspace(1,40,40);%voxels
test_1((ii-1)*40+1:(ii-1)*40+40,4) = [1,1,1,2,2,2,2,2,2,2,2,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6];% layers
end
T_1 = array2table(test_1,'VariableNames',{'lagtime','coeff','voxels','layers'});% cvt 2 table
label_1 = table2array(T_1(:,4));% Groups: labelling along layers
%% Original plots
g = figure; 
set(g, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['Original 3d scatter']);
h=scatter3(test_1(:,1),test_1(:,2),test_1(:,3),25,label_1,'filled');
%h.MarkerFaceColor = [0 0.5 0.5];
xlabel('lag time')
xlim([-10,10]);
ylabel('xc coeff')
zlabel('Voxels')

% statistics
g = figure; 
laglim = 30;
set(g, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['histogram']);
h_1 = histogram(test_1(:,1),'Normalization','pdf');
h_1.NumBins=50,h_1.FaceColor = [0 0.4470 0.7410];
hold on;
pf_1 = fitdist(test_1(:,1),'Normal');
x_1 = [-laglim:.1:laglim];% lag time interval
y_1 = normpdf(x_1,pf_1.mu,pf_1.sigma);
plot(x_1,y_1,'--','LineWidth',2);
xlabel('lag time');
ylabel('probabilities summation');
title(['\sigma:',num2str(pf_1.sigma),'   \mu:',num2str(pf_1.mu)]);

% original 2d scatter
g = figure; 
set(g, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['Original 2d scatter']);
markersize = 25;
hold on;
for ii = 1:6
scatter(T_1.lagtime(T_1.layers==ii),T_1.coeff(T_1.layers==ii),markersize,'o','filled')
end
xlinetitle = strcat('\mu = ',num2str(pf_1.mu));
xl = xline(pf_1.mu,'--m',xlinetitle);
xl.LineWidth = 2.5;
xl.LabelVerticalAlignment = 'top';
xl.LabelHorizontalAlignment = 'center';
xl.LabelOrientation = 'horizontal';
hold off
xlim([-laglim,laglim]);
xlabel('lag time');
ylabel('xc coeff');
grid on;
title('xc coeff against lag time of the animal');

% colormap of coeff to evaluate boundary selection
for ii = 1:3
g = figure; 
set(g, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['colormap coeff']);
testt = pks(:,:,ii);imagesc(testt);xlabel('trails');ylabel('Cortical voxels');colorbar
end
% tSNE plots
g = figure; 
set(g, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['tSNE 3d scatter']);
perpl = 30;
Y3 = tsne(test_1,'Algorithm','exact','Distance','euclidean','NumPCAComponents',3,'NumDimensions',3,'Perplexity',perpl,'Standardize',true);
scatter3(Y3(:,1),Y3(:,2),Y3(:,3),25,label_1,'filled');
% 2D tSNE:
g = figure; 
set(g, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['tSNE 2d scatter']);
gs = gscatter(Y3(:,1),Y3(:,2),label_1);
view(-93,14)
% another hist
g = figure; 
set(g, 'WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', ['histogram']);
h_1 = histogram(T_1.lagtime);
h_1.NumBins=50,h_1.FaceColor = [0 0.4470 0.7410];


% [hc_1,ct_1] = histcounts(T_1.lagtime,50);
% [hc_1,idx_1] = sort(hc_1,'descend');
% for ii=1:50
%     ct_1(ii)=ct_1(idx_1(ii));
% end
% ct_1 = unique(ct_1);
% bh_1 = bar(ct_1,hc_1(1:length(unique(ct_1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ABORTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PCA plots
% [coeff_1,score_1,latent_1,tsquared_1,explained_1] = pca(locs(:,:,3));
% h=scatter3(score_1(:,1),score_1(:,2),score_1(:,3),25,label_1,'filled');
% %h.MarkerFaceColor = [0 0.5 0.5];
% axis equal
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')
% % PCA colormap
% [coeff_1,score_1,latent_1] = pca(pks(:,:,3));
% Xcentered_1 = score_1*coeff_1'
% imagesc(Xcentered_1);colorbar