if savefigs
    close all
    clear fig_hand
    fig_number = 1;
end

%% Dependencies
% The following code relies on the Trees toolbox (https://www.treestoolbox.org/)
% as well as these dependencies which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/
%   distinguishable_colors 
%   BeeswarmPlot 
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
filter_ODI_logical = ismember(filter_Morph,filter_Morph_ODI);
filter_ODI_idx = find(filter_ODI_logical);

pca_based_elongation_measure = cellfun(@(x) x.elongation_ratio,{data(filter_Morph).morphology});
centerofmas_based_asymetry_measure = cellfun(@(x) x.asymetry_ratio,{data(filter_Morph).morphology});
sholl_max_crossings = cellfun(@(x) max(x.Sholl_analysis.sholl_analysis.sd),{data(filter_Morph).morphology});
DOi = cellfun(@(x) x.DOi,{data(filter_Morph).morphology});

ODI = [data(filter_Morph_ODI).ODI_AMPA_step_peak];
eye_pref = ODI>0;
input_type = ODI; input_type(ODI>0.8) = 1; input_type(ODI<-0.8) = 3; input_type(abs(ODI)<0.8) = 2;

clear all_radial_2d_densities all_radial_3d_densities all2D_radial_d_densities dendrite_length
Distance_Rdensity =  1.61:1.61:300;

for i = 1:length(data(filter_Morph))
    try
        all_radial_3d_densities(i,:) = data(filter_Morph(i)).radial_3D_dendrite_density(1,:);
        dendrite_length(i,:) = data(filter_Morph(i)).morphology.Sholl_analysis.total_length;
        dendrite_stem_num(i,:) = data(filter_Morph(i)).morphology.Sholl_analysis.num_tips;
        
    catch
        all_radial_3d_densities(i,:) = nan(1,length(Distance_Rdensity));
        dendrite_length(i,:) = nan;
        dendrite_stem_num(i,:) = nan;
    end
end

dendrite_length(find(sum(all_radial_3d_densities,2)==0)) = [];
dendrite_stem_num(find(sum(all_radial_3d_densities,2)==0)) = [];
all_radial_3d_densities(find(sum(all_radial_3d_densities,2)==0),:) = [];

for i = 1:size(all_radial_3d_densities,1)
    max_dend_length(i)=data(filter_Morph(i)).morphology.Sholl_analysis.max_euclidean_distance;
end

Rdensity_mean = nanmean(all_radial_3d_densities,1)';

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'Example_morpholgies', 'Position', [1,716,1025,271]);
fig_number = fig_number+1;
rad_ex = '190206SetupB0003';
rad_ex_idx = find([ismember({data(filter_Morph).patching_date},[rad_ex(1:6) 'LGN']) &...
    ismember({data(filter_Morph).Setup},rad_ex(7:12)) &...
    ismember({data(filter_Morph).cellname},rad_ex(13:16))]);
rad_ex_idx2 = find(filter_ODI_idx==rad_ex_idx);
ax(1) = subplot(1,3,1); %subplot(2,2,1.5);
scatter(data(filter_Morph(rad_ex_idx)).morphology.interp_morphology_transformed(:,1),...
    data(filter_Morph(rad_ex_idx)).morphology.interp_morphology_transformed(:,2),2,'k','filled'); hold on
scatter(0,0,50,'m','filled')
elongation_ratio=data(filter_Morph(rad_ex_idx)).morphology.elongation_ratio;
elongation_angle=data(filter_Morph(rad_ex_idx)).morphology.elongation_angle;
asymetry_ratio=data(filter_Morph(rad_ex_idx)).morphology.asymetry_ratio;
asymetry_angle=data(filter_Morph(rad_ex_idx)).morphology.asymetry_angle;
max_dend_length_temp=max(sqrt(data(filter_Morph(rad_ex_idx)).morphology.interp_morphology_transformed(:,1).^2+...
    data(filter_Morph(rad_ex_idx)).morphology.interp_morphology_transformed(:,1).^2));
[elong_cat1, elong_cat2] = pol2cart(elongation_angle,elongation_ratio*max_dend_length_temp);
[asym_cat1, aym_cat2] = pol2cart(asymetry_angle,asymetry_ratio*max_dend_length_temp);
plot([-elong_cat1 elong_cat1],[-elong_cat2 elong_cat2],'--m')
plot([0 asym_cat1],[0 aym_cat2],'m')
axis square off; box off;
title({'Y';['ODI: ' num2str(data(filter_Morph(rad_ex_idx)).ODI_AMPA_step_peak,2)];...
    ['DOi: ' num2str(data(filter_Morph(rad_ex_idx)).morphology.DOi,2)]})

elong_ex = '181219SetupA0004';
elong_ex_idx = find([ismember({data(filter_Morph).patching_date},[elong_ex(1:6) 'LGN']) &...
    ismember({data(filter_Morph).Setup},elong_ex(7:12)) &...
    ismember({data(filter_Morph).cellname},elong_ex(13:16))]);
elong_ex_idx2 = find(filter_ODI_idx==elong_ex_idx);
ax(2) = subplot(1,3,2);
scatter(data(filter_Morph(elong_ex_idx)).morphology.interp_morphology_transformed(:,1),...
    data(filter_Morph(elong_ex_idx)).morphology.interp_morphology_transformed(:,2),2,'k','filled'); hold on
scatter(0,0,50,'m','filled')
elongation_ratio=data(filter_Morph(elong_ex_idx)).morphology.elongation_ratio;
elongation_angle=data(filter_Morph(elong_ex_idx)).morphology.elongation_angle;
asymetry_ratio=data(filter_Morph(elong_ex_idx)).morphology.asymetry_ratio;
asymetry_angle=data(filter_Morph(elong_ex_idx)).morphology.asymetry_angle;
max_dend_length_temp=max(sqrt(data(filter_Morph(elong_ex_idx)).morphology.interp_morphology_transformed(:,1).^2+...
    data(filter_Morph(elong_ex_idx)).morphology.interp_morphology_transformed(:,1).^2));
[elong_cat1, elong_cat2] = pol2cart(elongation_angle,elongation_ratio*max_dend_length_temp);
[asym_cat1, aym_cat2] = pol2cart(asymetry_angle,asymetry_ratio*max_dend_length_temp);
plot([-elong_cat1 elong_cat1],[-elong_cat2 elong_cat2],'--m')
plot([0 asym_cat1],[0 aym_cat2],'m')
axis square off;box off;
title({'X';['ODI: ' num2str(data(filter_Morph(elong_ex_idx)).ODI_AMPA_step_peak,2)];...
    ['DOi: ' num2str(data(filter_Morph(elong_ex_idx)).morphology.DOi,2)]})

asym_ex = '181217SetupA0005';
asym_ex_idx = find([ismember({data(filter_Morph).patching_date},[asym_ex(1:6) 'LGN']) &...
    ismember({data(filter_Morph).Setup},asym_ex(7:12)) &...
    ismember({data(filter_Morph).cellname},asym_ex(13:16))]);
asym_ex_idx2 = find(filter_ODI_idx==asym_ex_idx);
ax(3) = subplot(1,3,3); 
scatter(data(filter_Morph(asym_ex_idx)).morphology.interp_morphology_transformed(:,1),...
    data(filter_Morph(asym_ex_idx)).morphology.interp_morphology_transformed(:,2),2,'k','filled'); hold on
scatter(0,0,50,'m','filled')
elongation_ratio=data(filter_Morph(asym_ex_idx)).morphology.elongation_ratio;
elongation_angle=data(filter_Morph(asym_ex_idx)).morphology.elongation_angle;
asymetry_ratio=data(filter_Morph(asym_ex_idx)).morphology.asymetry_ratio;
asymetry_angle=data(filter_Morph(asym_ex_idx)).morphology.asymetry_angle;
max_dend_length_temp=max(sqrt(data(filter_Morph(asym_ex_idx)).morphology.interp_morphology_transformed(:,1).^2+...
    data(filter_Morph(asym_ex_idx)).morphology.interp_morphology_transformed(:,1).^2));
[elong_cat1, elong_cat2] = pol2cart(elongation_angle,elongation_ratio*max_dend_length_temp);
[asym_cat1, aym_cat2] = pol2cart(asymetry_angle,asymetry_ratio*max_dend_length_temp);
plot([-elong_cat1 elong_cat1],[-elong_cat2 elong_cat2],'--m')
plot([0 asym_cat1],[0 aym_cat2],'m')
axis square off;box off;
title({'W';['ODI: ' num2str(data(filter_Morph(asym_ex_idx)).ODI_AMPA_step_peak,2)];...
    ['DOi: ' num2str(data(filter_Morph(asym_ex_idx)).morphology.DOi,2)]})
linkaxes(ax(:)); 
plot([max(xlim)-100 max(xlim)],[min(ylim) min(ylim)],'k','LineWidth',3)

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'Elongation_vs_asymmetry', 'Position', [1039, 717, 434, 268]); 
fig_number = fig_number+1;
scatter(pca_based_elongation_measure(:),centerofmas_based_asymetry_measure(:),'filled',...
    'MarkerFaceColor','k','MarkerEdgeColor','w'...
    ,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); hold on
% scatter(pca_based_elongation_measure(filter_ODI_idx(eye_pref)),...
%     centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref)),'b','filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); hold on
% scatter(pca_based_elongation_measure(filter_ODI_idx(~eye_pref)),...
%     centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref)),'r','filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); 
scatter(pca_based_elongation_measure([rad_ex_idx elong_ex_idx asym_ex_idx]),...
    centerofmas_based_asymetry_measure([rad_ex_idx elong_ex_idx asym_ex_idx]),[],'k','filled','MarkerEdgeColor','k');
xlabel(['Elongation']); ylabel(['Asymmetry'])
xlim([0 1]); ylim([0 1]); hold off
xticks([0:0.2:1]); yticks([0:0.2:1]);
axis square
[r,p] = corr(pca_based_elongation_measure(:),centerofmas_based_asymetry_measure(:));
title(['n=' num2str(length(pca_based_elongation_measure)), ', r=' num2str(r,2) ', p=' num2str(p,2)])

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'DOi', 'Position', [1,329,590,297]); 
fig_number = fig_number+1;
ax = subplot(1,2,1);
histogram(DOi,[0:0.1:1],'FaceColor',[0.5 0.5 0.5])
axis square
xlim([0 1]); xticks([0:0.2:1]); box off
xlabel('DOi'); ylabel('count')
originalSize = get(ax, 'Position');
set(ax, 'Position',originalSize-[0.02 0 0 0]);

ax = subplot(1,2,2);
scatter(pca_based_elongation_measure(:),centerofmas_based_asymetry_measure(:),[],DOi,'filled'); hold on
xlabel(['Elongation']); ylabel(['Asymmetry'])
xlim([0 1]); ylim([0 1]); xticks([0:0.2:1]); yticks([0:0.2:1]);
axis square
originalSize = get(ax, 'Position');
ac = colorbar;
ac.Title.String = 'DOi';
ac.TickDirection = 'out';
set(ax, 'Position',originalSize-[0.02 0 0 0]);

Conf_slice_align_angle = cellfun(@(x) wrapToPi(deg2rad(x)),{data(filter_Morph_ccf_alinged).Conf_to_CCF_rotation});
elongation_angle = cellfun(@(x) x.elongation_angle,{data(filter_Morph_ccf_alinged).morphology});
centermass_angle = cellfun(@(x) x.asymetry_angle,{data(filter_Morph_ccf_alinged).morphology});

elongation_angle = wrapToPi((elongation_angle+Conf_slice_align_angle).*2)./2;
centermass_angle = wrapToPi((centermass_angle+Conf_slice_align_angle));

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'asym dir and elong ori', 'Position', [1,41,934,295]); 
fig_number = fig_number+1;
clear ax gca
subplot(1,3,1); polarhistogram(elongation_angle,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(1) = gca; 
hold on; polarplot(repmat(circ_mean(elongation_angle'.*2)./2,1,2),(rlim),'--k','LineWidth',3)
title({'elong. ori'; ['mean: ' num2str(wrapTo180(rad2deg(circ_mean(elongation_angle'.*2)./2)),3) '°'];...
    ['Rayleigh test p: ' num2str(circ_rtest(elongation_angle),2)]});
subplot(1,3,2); polarhistogram(centermass_angle,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(2) = gca; 
hold on; polarplot(repmat(circ_mean(centermass_angle'),1,2),(rlim),'--k','LineWidth',3)
title({'asym. dir';['mean: ' num2str(wrapTo180(rad2deg(circ_mean(centermass_angle'))),3) '°'];...
    ['Rayleigh test p: ' num2str(circ_rtest(centermass_angle),2)]}); 
ax(1).ThetaLim=[-90 90]; ax(1).ThetaTick=(-90:90:90); ax(1).ThetaTickLabels={'-90°','0°','90°'};
ax(2).ThetaTick=(0:90:360); ax(2).ThetaTickLabels={'0°','90°','+/-180°','-90°'};
for i = 1:length(ax)
    ax(i).ThetaDir='clockwise';
    ax(i).RAxisLocation = 270;
end
subplot(1,3,3)
scatter(wrapToPi((elongation_angle-circ_mean(elongation_angle'))*2)/2,wrapToPi(centermass_angle-circ_mean(centermass_angle')))
[R,p]=circ_corrcc(elongation_angle,centermass_angle);
title(['corr. R=' num2str(R,2) ', p=' num2str(p,2)])
xlabel('elongation angle')
ylabel('asymmetry angle')
axis tight
xticks([-pi/2:pi/4:pi/2]); xticklabels({'-90','-45','0','45','90'})
yticks([-pi:pi/2:pi]); yticklabels({'-180','-90','0','90','180'})

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'IvsC_cell_morpho', 'Position', [260,49,1545,508]);
fig_number = fig_number+1;
subplot(2,5,1)
plotSpread({max_dend_length(filter_ODI_idx(eye_pref)),max_dend_length(filter_ODI_idx(~eye_pref))},...
    'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(max_dend_length(filter_ODI_idx(eye_pref))),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(max_dend_length(filter_ODI_idx(~eye_pref))),1,2),'k','LineWidth',3)
xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('max dend reach'); axis square
temp1 = max_dend_length(filter_ODI_idx(eye_pref));
temp2 = max_dend_length(filter_ODI_idx(~eye_pref));
[kpcontra,~]=kstest((temp1-mean(temp1))/std(temp1));
[kpipsi,~]=kstest((temp2-mean(temp2))/std(temp1));
if ~kpcontra & ~kpipsi
    [~,p,~,STATS]=ttest2(temp1,temp2,'tail','both','vartype','unequal');
    title({'ttest'; [['T = ' num2str(STATS.tstat(1),3)], [', p = ' num2str(p,3)], ['df = ' num2str(STATS.df(1),3)]]})
else
    [p,~,STATS]=ranksum(temp1,temp2);
    title({'Wilcoxon RS test'; [['W= ' num2str(STATS.ranksum)],[', p= ' num2str(p,3)]]})
end
ylim([0 inf])

subplot(2,5,2)
plotSpread({dendrite_length(filter_ODI_idx(eye_pref)),dendrite_length(filter_ODI_idx(~eye_pref))},...
    'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(dendrite_length(filter_ODI_idx(eye_pref))),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(dendrite_length(filter_ODI_idx(~eye_pref))),1,2),'k','LineWidth',3)
xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('total dend length'); axis square
temp1 = dendrite_length(filter_ODI_idx(eye_pref));
temp2 = dendrite_length(filter_ODI_idx(~eye_pref));
[kpcontra,~]=kstest((temp1-mean(temp1))/std(temp1));
[kpipsi,~]=kstest((temp2-mean(temp2))/std(temp1));
if ~kpcontra & ~kpipsi
    [~,p,~,STATS]=ttest2(temp1,temp2,'tail','both','vartype','unequal');
    title({'ttest'; [['T = ' num2str(STATS.tstat(1),3)], [', p = ' num2str(p,3)], ['df = ' num2str(STATS.df(1),3)]]})
else
    [p,~,STATS]=ranksum(temp1,temp2);
    title({'Wilcoxon RS test'; [['W= ' num2str(STATS.ranksum)],[', p= ' num2str(p,3)]]})
end

subplot(2,5,3)
plotSpread({sholl_max_crossings(filter_ODI_idx(eye_pref)),sholl_max_crossings(filter_ODI_idx(~eye_pref))},...
    'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(sholl_max_crossings(filter_ODI_idx(eye_pref))),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(sholl_max_crossings(filter_ODI_idx(~eye_pref))),1,2),'k','LineWidth',3)
xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('sholl max crossing'); axis square
temp1 = sholl_max_crossings(filter_ODI_idx(eye_pref));
temp2 = sholl_max_crossings(filter_ODI_idx(~eye_pref));
[kpcontra,~]=kstest((temp1-mean(temp1))/std(temp1));
[kpipsi,~]=kstest((temp2-mean(temp2))/std(temp1));
if ~kpcontra & ~kpipsi
    [~,p,~,STATS]=ttest2(temp1,temp2,'tail','both','vartype','unequal');
    title({'ttest'; [['T = ' num2str(STATS.tstat(1),3)], [', p = ' num2str(p,3)], ['df = ' num2str(STATS.df(1),3)]]})
else
    [p,~,STATS]=ranksum(temp1,temp2);
    title({'Wilcoxon RS test'; [['W= ' num2str(STATS.ranksum)],[', p= ' num2str(p,3)]]})
end
ylim([0 inf])

subplot(2,5,6); 
plotSpread({DOi(filter_ODI_idx(eye_pref)),DOi(filter_ODI_idx(~eye_pref))},...
    'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(DOi(filter_ODI_idx(eye_pref))),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(DOi(filter_ODI_idx(~eye_pref))),1,2),'k','LineWidth',3)
xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('DOi'); ylim([0 1]); axis square; 
temp1 = DOi(filter_ODI_idx(eye_pref));
temp2 = DOi(filter_ODI_idx(~eye_pref));
[kpcontra,~]=kstest((temp1-mean(temp1))/std(temp1));
[kpipsi,~]=kstest((temp2-mean(temp2))/std(temp1));
if ~kpcontra & ~kpipsi
    [~,p,~,STATS]=ttest2(temp1,temp2,'tail','both','vartype','unequal');
    title({'ttest'; [['T = ' num2str(STATS.tstat(1),3)], [', p = ' num2str(p,3)], ['df = ' num2str(STATS.df(1),3)]]})
else
    [p,~,STATS]=ranksum(temp1,temp2);
    title({'Wilcoxon RS test'; [['W= ' num2str(STATS.ranksum)],[', p= ' num2str(p,3)]]})
end
ylim([0 inf])

subplot(2,5,7)
plotSpread({centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref)),centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref))},...
    'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref))),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref))),1,2),'k','LineWidth',3)
xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('asymetry'); axis square
temp1 = centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref));
temp2 = centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref));
[kpcontra,~]=kstest((temp1-mean(temp1))/std(temp1));
[kpipsi,~]=kstest((temp2-mean(temp2))/std(temp1));
if ~kpcontra & ~kpipsi
    [~,p,~,STATS]=ttest2(temp1,temp2,'tail','both','vartype','unequal');
    title({'ttest'; [['T = ' num2str(STATS.tstat(1),3)], [', p = ' num2str(p,3)], ['df = ' num2str(STATS.df(1),3)]]})
else
    [p,~,STATS]=ranksum(temp1,temp2);
    title({'MWU'; [['U= ' num2str(STATS.ranksum)],[', p= ' num2str(p,3)]]})
end

subplot(2,5,8)
plotSpread({pca_based_elongation_measure(filter_ODI_idx(eye_pref)),pca_based_elongation_measure(filter_ODI_idx(~eye_pref))},...
    'categoryLabels',{'contra','ipsi'},'distributionColors',{'b','r'})
plot([1-0.2 1+0.2],repmat(mean(pca_based_elongation_measure(filter_ODI_idx(eye_pref))),1,2),'k','LineWidth',3); hold on
plot([2-0.2 2+0.2],repmat(mean(pca_based_elongation_measure(filter_ODI_idx(~eye_pref))),1,2),'k','LineWidth',3)
xticks([1 2]); xticklabels({'contra' 'ipsi'})
ylabel('elongation'); axis square
temp1 = pca_based_elongation_measure(filter_ODI_idx(eye_pref));
temp2 = pca_based_elongation_measure(filter_ODI_idx(~eye_pref));
[kpcontra,~]=kstest((temp1-mean(temp1))/std(temp1));
[kpipsi,~]=kstest((temp2-mean(temp2))/std(temp1));
if ~kpcontra & ~kpipsi
    [~,p,~,STATS]=ttest2(temp1,temp2,'tail','both','vartype','unequal');
    title({'ttest'; [['T = ' num2str(STATS.tstat(1),3)], [', p = ' num2str(p,3)], ['df = ' num2str(STATS.df(1),3)]]})
else
    [p,~,STATS]=ranksum(temp1,temp2);
    title({'Wilcoxon RS test'; [['W= ' num2str(STATS.ranksum,3)],[', p= ' num2str(p,3)]]})
end

subplot(2,5,[4,5,9,10])
scatter(pca_based_elongation_measure(filter_ODI_idx(eye_pref)),...
    centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref)),'b','filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); hold on
scatter(pca_based_elongation_measure(filter_ODI_idx(~eye_pref)),...
    centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref)),'r','filled','MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3); 
elong_median_C = median(pca_based_elongation_measure(filter_ODI_idx(eye_pref)));
elong_p1_C = prctile(pca_based_elongation_measure(filter_ODI_idx(eye_pref)),25); 
elong_p2_C = prctile(pca_based_elongation_measure(filter_ODI_idx(eye_pref)),75); 
asym_median_C = mean(centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref)));
asym_p1_C = prctile(centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref)),25); 
asym_p2_C = prctile(centerofmas_based_asymetry_measure(filter_ODI_idx(eye_pref)),75); 
elong_median_I = median(pca_based_elongation_measure(filter_ODI_idx(~eye_pref)));
elong_p1_I = prctile(pca_based_elongation_measure(filter_ODI_idx(~eye_pref)),25); 
elong_p2_I = prctile(pca_based_elongation_measure(filter_ODI_idx(~eye_pref)),75); 
asym_median_I = mean(centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref)));
asym_p1_I = prctile(centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref)),25); 
asym_p2_I = prctile(centerofmas_based_asymetry_measure(filter_ODI_idx(~eye_pref)),75);  
errorbar(elong_median_C,asym_median_C,...
    asym_p1_C-asym_median_C,asym_p2_C-asym_median_C,...
    elong_p1_C-elong_median_C,elong_p2_C-elong_median_C,'color','b','LineWidth',3)
errorbar(elong_median_I,asym_median_I,...
    asym_p1_I-asym_median_I,asym_p2_I-asym_median_I,...
    elong_p1_I-elong_median_I,elong_p2_I-elong_median_I,'color','r','LineWidth',3)
xlabel(['Elongation']); ylabel(['Asymmetry'])
xlim([0 1]); ylim([0 1]); hold off
xticks([0:0.2:1]); yticks([0:0.2:1]);
axis square
title(['contra=' num2str(length(filter_ODI_idx(eye_pref))) ', ipsi=' num2str(length(filter_ODI_idx(~eye_pref)))])

%% save figs
if savefigs
    mkdir([savefigs_location '\morphfigs\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[savefigs_location '\morphfigs\' fig_hand(i).Name],'svg')
        saveas(fig_hand(i),[savefigs_location '\morphfigs\' fig_hand(i).Name],'fig')
    end
end

