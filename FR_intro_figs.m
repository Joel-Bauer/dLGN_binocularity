if savefigs
    close all
    clear fig_hand
    fig_number = 1;
end

%% Dependencies
% The following code relies on HartigansDipTest.m and HartigansDipSignifTest.m 
% as well as these dependencies which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/ 
%   PlotConfInts
%   BeeswarmPlot
%   cbrewer
%	distributionPlot
%   suptitle

%% Morph radial histogram
filter_Morph = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', nan, ...
    'dLGN', 1, ...
    'A_ramp_signal_QC',nan,...
    'CS_internal',nan,...
    'Interneuron_morph',0,...
    'morph_QC',1));

clear all_radial_3d_densities max_dend_length
Distance_Rdensity =  1.61:1.61:300;
for i = 1:length(data(filter_Morph))
    try
        all_radial_3d_densities(i,:) = data(filter_Morph(i)).radial_3D_dendrite_density(1,:);        
    catch
        all_radial_3d_densities(i,:) = nan(1,length(Distance_Rdensity));
    end
end
for i = 1:size(all_radial_3d_densities,1)
    max_dend_length(i)=data(filter_Morph(i)).morphology.Sholl_analysis.max_euclidean_distance;
end

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'Radial_histogram', 'Position', [1,701,746,294]); 
fig_number = fig_number+1;
subplot(1,2,1)
for i = 1:size(all_radial_3d_densities,1)
    h = plot(Distance_Rdensity,all_radial_3d_densities(i,:),'color',[0.6 0.6 0.6]); hold on
    h.Color(4)=0.2;
end
all_radial_3d_densities_cuttoff = all_radial_3d_densities;
all_radial_3d_densities_cuttoff(find(all_radial_3d_densities(:,end)~=0),:) = [];
Rdensity_mean = nanmean(all_radial_3d_densities_cuttoff,1)';
plot(Distance_Rdensity,Rdensity_mean,'k')
xlim([0 max(max_dend_length)]); box off; axis square
xlabel('distance from soma (um)')
ylabel('Fraction of total dendritic arbour')

subplot(1,2,2)
cum_Rdense = cumsum(Rdensity_mean);
plot(Distance_Rdensity,cum_Rdense,'k'); hold on
idx_75 = find(abs(cum_Rdense-0.75)==min(abs(cum_Rdense-0.75)));
plot(Distance_Rdensity([1 idx_75]),cum_Rdense([idx_75 idx_75]),'--k');
plot(Distance_Rdensity([idx_75 idx_75]),cum_Rdense([1 idx_75]),'--k');
xlabel('distance from soma (um)'); box off; axis square
ylabel('Cumulative fraction of total dendritic arbour')
title({'75% of dendritic length is within'; [num2str(Distance_Rdensity(idx_75),2) 'um of the soma']})


for i = 1: length(Distance_Rdensity)
    if i == 1
        Distance_bin_area(i) = pi*4/3*(Distance_Rdensity(i)^3);
    else
        Distance_bin_area(i) = pi*4/3*(Distance_Rdensity(i)^3 - Distance_Rdensity(i-1)^3);
    end 
end

Rweightsum_mean = Rdensity_mean./Distance_bin_area'; 
Rweightsum_mean = Rweightsum_mean/sum(Rweightsum_mean);
Rweightsum_cumsum = cumsum(Rweightsum_mean);

figure
plot(Distance_Rdensity,Rweightsum_mean,'k')
xlim([0 max(max_dend_length)]); box off; axis square
xlabel('distance from soma (um)')
ylabel('Fraction of total dendritic arbour in 3D')

%% exampel method for Fr extraction 
example_cell = 22;
data_ex_cell = data(example_cell);

% mFR colorscheme
FR_colors = flip(cbrewer('div','RdYlGn',100));
% FR_colors = buildcmap('myc');

% top view
clear ratio_spher ratio_spher_side
ratio_spher(:,:,1) = nanmean(data_ex_cell.Fr_sphere_stacks.Stack_red,3);
ratio_spher(:,:,2) = nanmean(data_ex_cell.Fr_sphere_stacks.Stack_green,3);
ratio_spher(:,:,3) = zeros(size(ratio_spher,1),size(ratio_spher,2));
ratio_spher(isnan(ratio_spher)) = 0;

% side view
ratio_spher_side(:,:,1) = squeeze(nanmean(data_ex_cell.Fr_sphere_stacks.Stack_red,1))';
ratio_spher_side(:,:,2) = squeeze(nanmean(data_ex_cell.Fr_sphere_stacks.Stack_green,1))';
ratio_spher_side(:,:,3) = zeros(size(ratio_spher_side,1),size(ratio_spher_side,2));
ratio_spher_side(isnan(ratio_spher_side)) = 0;

Iratio_stack = (data_ex_cell.Fr_sphere_stacks.Stack_red - data_ex_cell.Fr_sphere_stacks.Stack_green)./(data_ex_cell.Fr_sphere_stacks.Stack_red + data_ex_cell.Fr_sphere_stacks.Stack_green);
morphology_transfomred_scalled = data_ex_cell.FlourRatio_data.morphology_transformed_scalled;

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'FR_calc_example', 'Position', [67 54 1773 914]);
fig_number = fig_number+1; clear ax
ax(1) = subplot(2,6,1); cla
h = imshow((ratio_spher)./255); hold on;
set(h,'AlphaData',mean(ratio_spher,3)>0)
scatter([morphology_transfomred_scalled(2,:)],...
    [morphology_transfomred_scalled(1,:)],5,'c','filled')
scatter(size(ratio_spher,1)/2,...
    size(ratio_spher,2)/2,40,'ob','filled');
daspect([4 4 4]); 
xticks([]); yticks([]); axis on
set(gca,'Color',[1 1 1])
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

ax(7) = subplot(2,6,7); cla
h = imshow((ratio_spher_side)./255); hold on; 
set(h,'AlphaData',mean(ratio_spher_side,3)>0)
scatter([morphology_transfomred_scalled(2,:)],...
    [morphology_transfomred_scalled(3,:)],5,'c','filled')
scatter(size(ratio_spher_side,2)/2, ...
    data_ex_cell.Confocal_cell_loc(3), 40,'ob','filled');
daspect([4 1.6 1.6]); 
xticks([]); yticks([]); axis on
set(gca,'Color',[1 1 1])
xlim(diff(xlim)/4*[1 -1]+(xlim)); 

Fr_binned = squeeze(nanmean(data_ex_cell.FlourRatio_data.Binned_ratio_image3D_all(:,:,:,1),3));
Morph_binned = full(squeeze(nanmean(data_ex_cell.FlourRatio_data.ThreeD_dendrite_density_rotated(:,:,:,1),3)));

ax(2) = subplot(2,6,2);
h=imagesc(Fr_binned); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
set(h,'AlphaData',~isnan(Fr_binned))
colormap(ax(2),FR_colors)
caxis([-0.7 0.7])
xticks([]); yticks([]); axis on
set(gca,'Color',[1 1 1])
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

ax(3) = subplot(2,6,3);
h=imagesc(Morph_binned); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
set(h,'AlphaData',~isnan(Morph_binned))
colormap(ax(3),flip(colormap(ax(3),'gray')))
xticks([]); yticks([]); axis on
set(gca,'Color',[1 1 1])
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

Fr_binned = squeeze(nanmean(data_ex_cell.FlourRatio_data.Binned_ratio_image3D_all(:,:,:,1),1));
Morph_binned = full(squeeze(nanmean(data_ex_cell.FlourRatio_data.ThreeD_dendrite_density_rotated(:,:,:,1),1)));
Morph_binned(Morph_binned==0) = nan;

ax(8) = subplot(2,6,8);
h=imagesc(Fr_binned'); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(2),data_ex_cell.FlourRatio_data.Binned_cell_pos(3),[],'b','filled')
set(h,'AlphaData',~isnan(Fr_binned'))
colormap(ax(8),FR_colors)
caxis([-0.7 0.7])
xticks([]); yticks([]); axis on
set(gca,'Color',[1 1 1])
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

ax(9) = subplot(2,6,9);
h=imagesc(Morph_binned'); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(2),data_ex_cell.FlourRatio_data.Binned_cell_pos(3),[],'b','filled')
set(h,'AlphaData',~isnan(Morph_binned'))
colormap(ax(9),flip(colormap(ax(9),'gray')))
xticks([]); yticks([]); axis on
set(gca,'Color',[1 1 1])
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

Fr_binned = data_ex_cell.FlourRatio_data.Binned_ratio_image3D_all(:,:,:,1);
Morph_binned = full(data_ex_cell.FlourRatio_data.ThreeD_dendrite_density_rotated(:,:,:,1));
Morph_binned(Morph_binned==0) = nan;
Fr_binned(isnan(Morph_binned)) = nan; % apply morph mask
Fr_binned_zproj =  squeeze(nanmean((Fr_binned.*Morph_binned),3))./squeeze(nanmean(Morph_binned,3));
Fr_binned_xproj =  squeeze(nanmean((Fr_binned.*Morph_binned),1))./squeeze(nanmean(Morph_binned,1));
Morph_binned_zproj = squeeze(nanmean(Morph_binned,3));
Morph_binned_xproj = squeeze(nanmean(Morph_binned,1));

ax(4) = subplot(2,6,4); cla
h = imagesc(Fr_binned_zproj); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
set(h,'AlphaData',...
    (Morph_binned_zproj./max(Morph_binned_zproj(:)).*(1+0.2)+0.2)...
    .*~isnan(Morph_binned_zproj)...
    .*~isnan(Fr_binned_zproj));
caxis([-0.7 0.7])
xticks([]); yticks([]);
colormap(ax(4),FR_colors)
pbaspect([1 1 1]); set(gca,'TickDir','out'); %colorbar
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

ax(10) = subplot(2,6,10); cla
h = imagesc(Fr_binned_xproj'); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(2),data_ex_cell.FlourRatio_data.Binned_cell_pos(3),[],'b','filled')
set(h,'AlphaData',...
    (Morph_binned_xproj'./max(Morph_binned_xproj(:)).*(1-0.2)+0.2)...
    .*~isnan(Morph_binned_xproj')...
    .*~isnan(Fr_binned_xproj'));
% set(gca,'Color',[0.7 0.7 0.7])
caxis([-0.7 0.7])
xticks([]); yticks([]);
colormap(ax(10),FR_colors)
pbaspect([1 1 1]); set(gca,'TickDir','out'); %colorbar
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

% example cell with spherical mask

Fr_binned = data_ex_cell.FlourRatio_data.Binned_ratio_image3D_all(:,:,:,1);
Morph_binned = data_ex_cell.FlourRatio_data.MeanRadial_3d_mask(:,:,:);
Fr_binned(isnan(Morph_binned)) = nan; % apply morph mask
Fr_binned_zproj =  squeeze(nanmean((Fr_binned.*Morph_binned),3))./squeeze(nanmean(Morph_binned,3));
Fr_binned_xproj =  squeeze(nanmean((Fr_binned.*Morph_binned),1))./squeeze(nanmean(Morph_binned,1));
Morph_binned(Morph_binned==0) = nan;
% Morph_binned(isnan(Fr_binned)) = nan;
Morph_binned_zproj = squeeze(nanmean(Morph_binned,3));
Morph_binned_xproj = squeeze(nanmean(Morph_binned,1));

ax(5) = subplot(2,6,5);
h=imagesc(Morph_binned_zproj); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
set(h,'AlphaData',~isnan(Morph_binned_zproj))
colormap(ax(5),flip(colormap(ax(5),'gray')))
xticks([]); yticks([]); axis on
set(gca,'Color',[1 1 1])
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

ax(6) = subplot(2,6,6); cla
h = imagesc(Fr_binned_zproj); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
set(h,'AlphaData',...
    (Morph_binned_zproj./max(Morph_binned_zproj(:)).*(1-0)+0)...
    .*~isnan(Morph_binned_zproj)...
    .*~isnan(Fr_binned_zproj));
% set(gca,'Color',[0.7 0.7 0.7])
caxis([-0.7 0.7])
colormap(ax(6),FR_colors)
xticks([]); yticks([]); axis on
set(gca,'TickDir','out'); %colorbar
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

ax(11) = subplot(2,6,11); cla
h=imagesc(Morph_binned_xproj'); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(2),data_ex_cell.FlourRatio_data.Binned_cell_pos(3),[],'b','filled')
set(h,'AlphaData',~isnan(Morph_binned_xproj'))
colormap(ax(11),flip(colormap(ax(11),'gray')))
xticks([]); yticks([]); axis on
set(gca,'TickDir','out'); %colorbar
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

ax(12) = subplot(2,6,12); cla
h = imagesc(Fr_binned_xproj'); hold on; axis image 
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(2),data_ex_cell.FlourRatio_data.Binned_cell_pos(3),[],'b','filled')
set(h,'AlphaData',...
    (Morph_binned_xproj'./max(Morph_binned_xproj(:)).*(1-0)+0)...
    .*~isnan(Morph_binned_xproj')...
    .*~isnan(Fr_binned_xproj'));
% set(gca,'Color',[0.7 0.7 0.7])
caxis([-0.7 0.7])
colormap(ax(12),FR_colors)
xticks([]); yticks([]); axis on
set(gca,'TickDir','out'); %colorbar
xlim(diff(xlim)/4*[1 -1]+(xlim)); ylim(diff(ylim)/4*[1 -1]+(ylim))

for i = 2:length(ax)
    ax(i).Position(3)=ax(1).Position(3);
end

%% general ipsi contra dom eye Fratios
filter_nonMD_Morph = find(QC_cell_filter(data, your_data_dir,...
    'transduction_QC', 1, 'dLGN', 1, ...
    'A_ramp_signal_QC', 1, 'N_ramp_signal_QC', nan, ...
    'Conf_QC', 1, 'TwoPTrans_QC', 1, 'CS_internal',1,'Interneuron_morph',0,...
    'morph_QC',1));
perc_morph_cut_norot = cell2mat(cellfun(@(x) x.Binned_ratio_images3D_percent_dendrites_excluded(1,1), {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false));
filter_nonMD_Morph(perc_morph_cut_norot>30) = [];

ODI_Morph = [data((filter_nonMD_Morph)).ODI_AMPA_step_peak];
eye_pref_Morph = ODI_Morph>0;

filter_nonMD_nonMorph = find(QC_cell_filter(data, your_data_dir,...
    'transduction_QC', 1, 'dLGN', 1, ...
    'A_ramp_signal_QC', 1, 'N_ramp_signal_QC', nan, ...
    'Conf_QC', 1, 'TwoPTrans_QC', nan, 'CS_internal',1,'Interneuron_morph',0,...
    'morph_QC',nan));
ODI_nonMorph = [data((filter_nonMD_nonMorph)).ODI_AMPA_step_peak];
eye_pref_nonMorph = ODI_nonMorph>0;

temp = cellfun(@(x) x.MorphDensity3D_adjusted_Fratio, {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false);
ThreeD_MorphDens = cat(1,temp{:});
temp = cellfun(@(x) x.MeanRadial3D_adjusted_Fratio, {data((filter_nonMD_nonMorph)).FlourRatio_data},'UniformOutput',false);
Radial_MorphDens = cat(1,temp{:});

radii = data(filter_nonMD_Morph(1)).FlourRatio_data.Sphere_radii;
sphere_size_used = 300;

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'mFR_performance', 'Position', [457, 464, 506, 530]);
fig_number = fig_number+1;
subplot(2,1,1)
contra_data = ThreeD_MorphDens(eye_pref_Morph,find(radii==sphere_size_used));
ipsi_data = ThreeD_MorphDens(~eye_pref_Morph,find(radii==sphere_size_used));
plotSpread({contra_data, ipsi_data}, 'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(contra_data),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(ipsi_data),1,2),'k','LineWidth',3)
set(gca,'TickDir','out'); xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('mFR'); axis square; ylim([-1 1])
[~,p,~,STATS]=ttest2(contra_data,ipsi_data,'tail','both','Vartype','unequal');
[~,~,shuff_output] = eyepref_decoder_boot(...
    [ThreeD_MorphDens(:,find(radii==sphere_size_used))],...
    ODI_Morph(:)',{'mFR'},0,0);
dec_ac_Morph = shuff_output.dec_acc(1);
title({['AMPA dom. p = ' num2str(p,2) ', t = ' num2str(STATS.tstat,2)];...
    ['n_c: ' num2str(sum(eye_pref_Morph),2) ', n_i: ' num2str(sum(~eye_pref_Morph),2) ]})

subplot(2,1,2)
histogram(shuff_output.dec_acc_shuf,30,'Normalization','probability'); hold on
percentile_grad = 95;
plot(repmat(prctile(shuff_output.dec_acc_shuf,percentile_grad),1,2),ylim,'--r')
plot(repmat(dec_ac_Morph,1,2),ylim,'--k')
xlabel('decoding accuracy'); xlim([0.4 1])
xticks([0.4:0.2:1])
ylabel('frequency')
box off; axis square
title({['decoding accuracy, black: ' num2str(dec_ac_Morph,2) ];...
    ['red: ' num2str(prctile(shuff_output.dec_acc_shuf,percentile_grad),2)]})

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'rFR_performance', 'Position', [966, 464, 506, 530]);
fig_number = fig_number+1;
subplot(2,1,1)
contra_data = Radial_MorphDens(eye_pref_nonMorph,find(radii==sphere_size_used));
ipsi_data = Radial_MorphDens(~eye_pref_nonMorph,find(radii==sphere_size_used));
plotSpread({contra_data, ipsi_data}, 'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(contra_data),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(ipsi_data),1,2),'k','LineWidth',3)
set(gca,'TickDir','out'); xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('rFR'); axis square; ylim([-1 1])
[~,p,~,STATS]=ttest2(contra_data,ipsi_data,'tail','both','Vartype','unequal');
[~,~,shuff_output] = eyepref_decoder_boot(...
    [Radial_MorphDens(:,find(radii==sphere_size_used))],...
    ODI_nonMorph(:)',{'rFR'},0,0);
dec_ac_nonMorph = shuff_output.dec_acc(1);
title({['AMPA dom. p = ' num2str(p,2) ', t = ' num2str(STATS.tstat,2)];...
    ['n_c: ' num2str(sum(eye_pref_nonMorph),2) ', n_i: ' num2str(sum(~eye_pref_nonMorph),2) ]})

subplot(2,1,2)
histogram(shuff_output.dec_acc_shuf,30,'Normalization','probability'); hold on
percentile_grad = 95;
plot(repmat(prctile(shuff_output.dec_acc_shuf,percentile_grad),1,2),ylim,'--r')
plot(repmat(dec_ac_nonMorph,1,2),ylim,'--k')
xlabel('decoding accuracy'); xlim([0.4 1])
xticks([0.4:0.2:1])
ylabel('frequency')
box off; axis square
title({['decoding accuracy, black: ' num2str(dec_ac_nonMorph,2) ];...
    ['red: ' num2str(prctile(shuff_output.dec_acc_shuf,percentile_grad),2)]})
%% direct comparison of mFR vs rFR
filter_nonMD_Morph = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', 1, 'dLGN', 1, ...
    'A_ramp_signal_QC', 1, 'N_ramp_signal_QC', nan, ...
    'Conf_QC', 1, 'TwoPTrans_QC', 1, 'CS_internal',1,'Interneuron_morph',0,...
    'morph_QC',1,'CCF_QC', 1));
perc_morph_cut_norot = cell2mat(cellfun(@(x) x.Binned_ratio_images3D_percent_dendrites_excluded(1,1), {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false));
perc_morph_cut_allrot_min = cell2mat(cellfun(@(x) min(x.Binned_ratio_images3D_percent_dendrites_excluded(:,1)), {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false));
perc_morph_cut_allrot_max = cell2mat(cellfun(@(x) max(x.Binned_ratio_images3D_percent_dendrites_excluded(:,1)), {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false));
perc_morph_cut_allrot_dif = perc_morph_cut_allrot_max-perc_morph_cut_allrot_min;
morph_asymetry = cell2mat(cellfun(@(x) x.asymetry_ratio, {data(filter_nonMD_Morph).morphology},'UniformOutput',false));

% 2 thresholds
filter_nonMD_Morph([perc_morph_cut_norot>30] | [perc_morph_cut_allrot_dif>10]) = [];
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'exclusion_threshold', 'Position', [1, 291, 453, 325]);
fig_number = fig_number+1;
scatter(perc_morph_cut_norot([perc_morph_cut_norot>30] | [perc_morph_cut_allrot_dif>10]),...
    perc_morph_cut_allrot_dif([perc_morph_cut_norot>30] | [perc_morph_cut_allrot_dif>10]),[],'b','filled',...
    'MarkerFaceAlpha',0.4); hold on
scatter(perc_morph_cut_norot(~[perc_morph_cut_norot>30] & ~[perc_morph_cut_allrot_dif>10]),...
    perc_morph_cut_allrot_dif(~[perc_morph_cut_norot>30] & ~[perc_morph_cut_allrot_dif>10]),[],'b','filled');
plot([min(xlim) max(xlim)],[10 10],'--k')
plot([30 30],[min(ylim) max(ylim)],'--k')
xlabel('% dend. cut at oritinal orientation')
ylabel({'max(delta(%)) dend. ' 'cut with morph rotation'})
axis square
title({['thresholds for exlcuding morphologies'] ['based on dendrites cut off']...
    ['included n: ' num2str(length(filter_nonMD_Morph),2)]})

ODI = [data((filter_nonMD_Morph)).ODI_AMPA_step_peak];
eye_pref = ODI>0;

Morph_rot_angles = data_ex_cell.FlourRatio_data.Morph_rot_angles;

% stats for delta Ori
shuffle_n = 10000;
boot_n = 1000;

temp = cellfun(@(x) x.MorphDensity3D_rotated_adjusted_Fratio(:,1), {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false);
ThreeD_MorphDens_rot = cat(2,temp{:})';
temp = cellfun(@(x) x.MeanRadial3D_adjusted_Fratio(:,1), {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false);
Radial_MorphDens_rot = cat(2,temp{:})';

[Dprime_p_value_matrix,~,shuff_output] = eyepref_decoder_boot(...
    [ThreeD_MorphDens_rot(:,1), ThreeD_MorphDens_rot(:,find(Morph_rot_angles==180)), Radial_MorphDens_rot],...
    ODI(:)',{'mFR' 'mFR180' 'rFR'},0,0);

mFR_dp = shuff_output.Dprime(1);
mFR_mean = mean(shuff_output.Dprime_bootstrp(1,:));
mFR_std = std(shuff_output.Dprime_bootstrp(1,:));

mFR180_dp = shuff_output.Dprime(2);
mFR180_mean = mean(shuff_output.Dprime_bootstrp(2,:));
mFR180_std = std(shuff_output.Dprime_bootstrp(2,:));

rFR_dp = shuff_output.Dprime(3);
rFR_mean = mean(shuff_output.Dprime_bootstrp(3,:));
rFR_std = std(shuff_output.Dprime_bootstrp(3,:));

mFR_vs_mFRrot_pval = Dprime_p_value_matrix(1,2);
mFR_vs_rFR_pval = Dprime_p_value_matrix(1,3);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'mFRvsrFR_dprime', 'Position', [1474, 698, 410, 295]);
fig_number = fig_number+1;
PlotConfInts( [1 2 3], [mFR_mean  rFR_mean mFR180_mean],...
    [mFR_std rFR_std mFR180_std], [0.8 0.8 0.8], 0, 'stderr', 'patch' ); hold on
plot([1 2 3], [mFR_mean rFR_mean mFR180_mean],'k')
scatter(1,mFR_dp,'^r','filled'); 
scatter(2,rFR_dp,'^r','filled');
scatter(3,mFR180_dp,'^r','filled');
ylim([0 2])
xlim([0.5 3.5])
xticks(1:3);
set(gca,'TickDir','out')
ylabel('dprime');
legend('off')
xticklabels({'mFR'  'rFR' 'mFR rot 180°'}); xtickangle(45)
axis square
box off
title({['mFR>mFR180° pval: ' num2str(mFR_vs_mFRrot_pval)];['mFR>rFR pval: ' num2str(mFR_vs_rFR_pval)]});

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'mFRvsmFRrotonly_dprime', 'Position', [1474, 698, 410, 295]);
fig_number = fig_number+1;
PlotConfInts( [1 2], [mFR_mean   mFR180_mean],...
    [mFR_std  mFR180_std], [0.8 0.8 0.8], 0, 'stderr', 'patch' ); hold on
plot([1 2], [mFR_mean  mFR180_mean],'k')
scatter(1,mFR_dp,'^r','filled'); 
scatter(2,mFR180_dp,'^r','filled');
ylim([0 2])
xlim([0.5 2.5])
xticks(1:2);
set(gca,'TickDir','out')
ylabel('dprime');
legend('off')
xticklabels({'mFR' 'mFR rot 180°'}); xtickangle(45)
axis square
box off
title({['mFR>mFR180° pval: ' num2str(mFR_vs_mFRrot_pval)]});

%% save figs
if savefigs
    mkdir([savefigs_location '\FRfigs\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[savefigs_location '\FRfigs\' fig_hand(i).Name],'svg')
        saveas(fig_hand(i),[savefigs_location '\FRfigs\' fig_hand(i).Name],'fig')
    end
end