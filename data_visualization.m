% data import
test_1(:,1) = locs(:,1,3);%lag times
test_1(:,2) = pks(:,1,3);%coeff
test_1(:,3) = linspace(1,40,40);%voxels
test_1(:,4) = [1,1,1,2,2,2,2,2,2,2,2,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6];% layers

T_1 = array2table(test_1,'VariableNames',{'lagtime','coeff','voxels','layers'});% cvt 2 table
label_1 = table2array(T_1(:,4));% labelling along layers
% Original plots
h=scatter3(test_1(:,1),test_1(:,2),test_1(:,3),25,label_1,'filled');
%h.MarkerFaceColor = [0 0.5 0.5];
xlabel('lag time')
ylabel('xc coeff')
zlabel('Voxels')

% PCA plots
[coeff_1,score_1,latent_1,tsquared_1,explained_1] = pca(locs(:,:,3));
h=scatter3(score_1(:,1),score_1(:,2),score_1(:,3),25,label_1,'filled');
%h.MarkerFaceColor = [0 0.5 0.5];
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
% PCA colormap
[coeff_1,score_1,latent_1] = pca(pks(:,:,3));
Xcentered_1 = score_1*coeff_1'
imagesc(Xcentered_1);colorbar

% tSNE plots
Y3 = tsne(test_1,'Algorithm','exact','Distance','euclidean','NumPCAComponents',3,'NumDimensions',3,'Perplexity',3,'Standardize',true);
figure
scatter3(Y3(:,1),Y3(:,2),Y3(:,3),25,label_1,'filled');
view(-93,14)
