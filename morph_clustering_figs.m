if savefigs
    close all
    clear fig_hand
    fig_number = 1;
end

%% Dependencies
% The following code relies on HartigansDipTest.m and
% HartigansDipSignifTest.m (http://www.nicprice.net/diptest/)
% as well as these dependencies which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/ 
%   cbrewer

%% morph 'types'. related
filter_Morph = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', nan, ...
    'dLGN', 1, ...
    'A_ramp_signal_QC',nan,...
    'CS_internal',nan,...
    'Interneuron_morph',0,...
    'morph_QC',1));
filter_Morph_ODI = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', 1, ...
    'dLGN', 1, ...
    'A_ramp_signal_QC',1,...
    'CS_internal',1,...
    'Interneuron_morph',0,...
    'morph_QC',1));
filter_Morph_ccf_alinged = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', nan, ...
    'dLGN', 1, ...
    'A_ramp_signal_QC',nan,...
    'CS_internal',nan,...
    'Interneuron_morph',0,...
    'Conf_QC', 1, ... %
    'TwoPTrans_QC', 1,... %
    'ccf', 1,... %
    'morph_QC',1));

%% Read out parameters
pca_based_elongation_measure = cellfun(@(x) x.elongation_ratio,{data(filter_Morph).morphology},'UniformOutput',0);
centerofmas_based_asymetry_measure = cellfun(@(x) x.asymetry_ratio,{data(filter_Morph).morphology},'UniformOutput',0);
for i = length(pca_based_elongation_measure)
   if isempty(pca_based_elongation_measure{i})
       pca_based_elongation_measure{i} = nan;
   end
   if isempty(centerofmas_based_asymetry_measure{i})
       centerofmas_based_asymetry_measure{i} = nan;
   end
end
pca_based_elongation_measure = cat(1,pca_based_elongation_measure{:});
centerofmas_based_asymetry_measure = cat(1,centerofmas_based_asymetry_measure{:});

sholl_max_crossings = cellfun(@(x) max(x.Sholl_analysis.sholl_analysis.sd),{data(filter_Morph).morphology});
DOi = cellfun(@(x) x.DOi,{data(filter_Morph).morphology});
clear max_dend_length total_l
for i = 1:size(filter_Morph,2)
    max_dend_length(i)=data(filter_Morph(i)).morphology.Sholl_analysis.max_euclidean_distance;
end
for i=1:length(filter_Morph)
     total_l(i)=[data(filter_Morph(i)).morphology.Sholl_analysis.total_length];
end

%% Combine parameters
morph_com=[max_dend_length' total_l' sholl_max_crossings' centerofmas_based_asymetry_measure pca_based_elongation_measure DOi'];
pred_lab={'Max Dendrite Length','Total Length','Max Sholl crossing','Asymetry','Elongation','DOi'};
nboot = 500;


%% dip test on morph var
% get the number of dimensions
data_dim = size(morph_com,2);

% allocate memory for the results
dip_results = zeros(data_dim,2);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'morph_variables_diptest', 'Position', [67,120,623,848]);
fig_number = fig_number+1;
% check for all PCs
for i = 1:data_dim
    [dip_results(i,1), dip_results(i,2)] = HartigansDipSignifTest(morph_com(:,i), nboot);
    
    % plot the result
%     ax(i) = subplot(ceil(sqrt(data_dim)),round(sqrt(data_dim)),i);
    ax(i) = subplot(round(sqrt(data_dim)),ceil(sqrt(data_dim)),i);
    if diff([min(morph_com(:,i)), max(morph_com(:,i))]) < 1 
        h1= histogram(morph_com(:,i),[0:0.1:1]);
        xlim([0 1])
    else
        h1= histogram(morph_com(:,i));
    end
    ylims(i,:) = ylim;
   
    box off;
    title(strjoin({'Dip test:',num2str(dip_results(i,1),2),'p val:',num2str(dip_results(i,2),2)},' '));
    ylabel('Count');
    xlabel(pred_lab{i})
end

for i = 1:data_dim
    ax(i).YLim=[0 round(max(ylims(:,2)),-1)]; 
    ax(i).YTick=[0:round(max(ylims(:,2)),-1)/4:round(max(ylims(:,2)),-1)];
end

%% morph var corr
[SpearR, corrp] = corr(morph_com,'type','Spearman');
bonf_corr_pthresh=0.05/(factorial(size(morph_com,2))/(factorial(size(morph_com,2)-2)*factorial(2)));

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'morph_variables_corr', 'Position', [699,176,1170,796]);
fig_number = fig_number+1;
a1 = subplot(2,3,[1,2,4,5]);
ax = imagesc(SpearR);
caxis([-1 1])
axis square
colormap(a1,cbrewer('div','RdBu',100));
% colormap('cool');
ax.Parent.XTickLabel = pred_lab;
ax.Parent.XAxisLocation = 'top';
ax.Parent.YTickLabel = pred_lab;
ax.Parent.XTick = 1:length(pred_lab);
ax.Parent.YTick = 1:length(pred_lab);
ax.Parent.XTickLabelRotation = -45;
ax.Parent.YTickLabelRotation = 0;
cb = colorbar;
cb.Title.String = 'Corr';

a2 = subplot(2,3,3);
imagesc(corrp<bonf_corr_pthresh)
colormap(a2,'gray'); caxis([0.5 0.6])
title({'sig corr 95% conf level'; 'with bonf corr'})
axis square
%% PCA
[coeff,score,latent,~,explained,mu] = pca([zscore(morph_com)]);

[pc_diptest(1,1), pc_diptest(1,2)]=HartigansDipSignifTest(score(:,1), nboot);
[pc_diptest(2,1), pc_diptest(2,2)]=HartigansDipSignifTest(score(:,2), nboot);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'PCA_dip_test', 'Position', [1150,373,546,598]);
fig_number = fig_number+1; 
h = scatterhist(score(:,1),score(:,2),'Kernel','off','Location','NorthEast',...
    'Direction','out','Color',[0.5 0.5 0.5],'LineStyle',{'-'},...
    'LineWidth',[2],'Marker','o','MarkerSize',[2]); hold on
plot([0 0],ylim,'k--')
plot(xlim,[0 0],'k--')
xlabel(['PC1, ' num2str(explained(1),2) '% exp. var.']); 
ylabel(['PC2, ' num2str(explained(2),2) '% exp. var.']); 
axis square
h(1).XAxisLocation= 'bottom' ;
h(1).YAxisLocation= 'left' ;
h(2).Visible = 'on';
h(2).Box = 'off';
h(2).Title.String = ['Dip test:',num2str(pc_diptest(1,1),2),', p val:',num2str(pc_diptest(1,2),2)];
h(2).YLim = [0 inf];
h(3).Visible = 'on';
h(3).Box = 'off';
h(3).Title.String = ['Dip test:',num2str(pc_diptest(2,1),2),', p val:',num2str(pc_diptest(2,2),2)];
h(3).YLim = [0 inf];

%% save figs
if savefigs
    mkdir([savefigs_location '\morph_antiClustering\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[savefigs_location '\morph_antiClustering\' fig_hand(i).Name],'svg')
        saveas(fig_hand(i),[savefigs_location '\morph_antiClustering\' fig_hand(i).Name],'fig')
    end
end