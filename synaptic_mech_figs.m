if savefigs
    close all
    clear fig_hand
    fig_number = 1;
end

%% Dependencies
% The following code relies on HartigansDipTest.m and HartigansDipSignifTest.m 
% as well as these dependencies which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/ 
%   distinguishable_colors
%   BeeswarmPlot

%% cells with similar Fr values can have oposite ODIs
% AMPA
filter_nonMD_Morph = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', 1, 'dLGN', 1, ...
    'A_ramp_signal_QC', 1, 'N_ramp_signal_QC', 1, ...
    'Conf_QC', 1, 'TwoPTrans_QC', 1, 'CS_internal',1,'Interneuron_morph',0,...
    'morph_QC',1));
perc_morph_cut_noro  t = cell2mat(cellfun(@(x) x.Binned_ratio_images3D_percent_dendrites_excluded(1,1), {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false));
filter_nonMD_Morph(perc_morph_cut_norot>30) = [];

example_cell = 22;
data_ex_cell = data(example_cell);

radii = data(filter_nonMD_Morph(1)).FlourRatio_data.Sphere_radii(1:end);
ODI_A = [data((filter_nonMD_Morph)).ODI_AMPA_step_peak];
ODI_N = [data((filter_nonMD_Morph)).ODI_NMDA_step_peak];

eye_pref_A = ODI_A>0;
eye_pref_N = ODI_N>0;

temp = cellfun(@(x) x.MorphDensity3D_adjusted_Fratio, {data((filter_nonMD_Morph)).FlourRatio_data},'UniformOutput',false);
ThreeD_MorphDens = cat(1,temp{:});

sphere_size_used = 300;
bins = 10;
binsOD  = -1:2/bins:1; % bin edges
% set(gca,'FontSize',13)

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'mFRvsODI_scatterplots', 'Position', [1   511   414   476]); 
fig_number = fig_number+1;

Fr_contra_histcount = histcounts(ThreeD_MorphDens(eye_pref_A,find(radii==sphere_size_used)),binsOD,'Normalization','count');
Fr_ipsi_histcount = histcounts(ThreeD_MorphDens(~eye_pref_A,find(radii==sphere_size_used)),binsOD,'Normalization','count');
CtI_ratio_over_Fr = (Fr_contra_histcount-Fr_ipsi_histcount)./(Fr_contra_histcount+Fr_ipsi_histcount);

subplot(2,1,1); 
scatter(ThreeD_MorphDens(eye_pref_A,find(radii==sphere_size_used)), ODI_A(eye_pref_A)',[],'b','filled','MarkerFaceAlpha',0.7); hold on
scatter(ThreeD_MorphDens(~eye_pref_A,find(radii==sphere_size_used)), ODI_A(~eye_pref_A)',[],'r','filled','MarkerFaceAlpha',0.7);
xlabel('axon fluor ratio'); xlim([-1 1])
ylabel({'AMPA based ODI'}); ylim([-1 1])
xticks(round(binsOD,2))
set(gca,'TickDir','out'); box off

Fr_contra_histcount = histcounts(ThreeD_MorphDens(eye_pref_N,find(radii==sphere_size_used)),binsOD,'Normalization','count');
Fr_ipsi_histcount = histcounts(ThreeD_MorphDens(~eye_pref_N,find(radii==sphere_size_used)),binsOD,'Normalization','count');
CtI_ratio_over_Fr = (Fr_contra_histcount-Fr_ipsi_histcount)./(Fr_contra_histcount+Fr_ipsi_histcount);

subplot(2,1,2);
scatter(ThreeD_MorphDens(eye_pref_N ,find(radii==sphere_size_used)), ODI_N(eye_pref_N )',[],'b','filled','MarkerFaceAlpha',0.7); hold on
scatter(ThreeD_MorphDens(~eye_pref_N ,find(radii==sphere_size_used)), ODI_N(~eye_pref_N )',[],'r','filled','MarkerFaceAlpha',0.7);
xlabel('axon fluor ratio'); xlim([-1 1])
ylabel({'NMDA based ODI'}); ylim([-1 1])
xticks(round(binsOD,2))
set(gca,'TickDir','out'); box off

%% AMPA to NMDA ratio
filter_nonMD_Morph = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', 1, 'dLGN', 1, ... 
    'A_ramp_signal_QC', 1, 'N_ramp_signal_QC', 1, 'AtoN_Rschange_QC',1,...
    'Conf_QC', nan, 'TwoPTrans_QC', nan, 'CS_internal', 1,'Interneuron_morph',0,...
    'morph_QC',nan));

ignore_missing_Nresp_cells = 1;
clear dom_AN_ratio nondom_AN_ratio dom_AN nondom_AN ODI_A ODI_N
for i = 1:length(filter_nonMD_Morph)
    if data(filter_nonMD_Morph(i)).brain_contra_ipsi == (data(filter_nonMD_Morph(i)).ODI_AMPA_step_peak>0)
        dom_AN_ratio(i) = data(filter_nonMD_Morph(i)).AMPA_NMDA_r_red;
        nondom_AN_ratio(i) = data(filter_nonMD_Morph(i)).AMPA_NMDA_r_blue;
    else
        dom_AN_ratio(i) = data(filter_nonMD_Morph(i)).AMPA_NMDA_r_blue;
        nondom_AN_ratio(i) = data(filter_nonMD_Morph(i)).AMPA_NMDA_r_red;
    end
    
    if isinf(dom_AN_ratio(i))
        if ignore_missing_Nresp_cells
            dom_AN(i) = nan;
        else
            dom_AN(i) = 1;
        end
    else
        dom_AN(i) = (dom_AN_ratio(i)-1)/(dom_AN_ratio(i)+1);
    end
    if isinf(nondom_AN_ratio(i))
        if ignore_missing_Nresp_cells
            nondom_AN(i) = nan;
        else
            nondom_AN(i) = 1;    
        end
    else
        nondom_AN(i) = (nondom_AN_ratio(i)-1)/(nondom_AN_ratio(i)+1);
    end
    ODI_A(i)=data(filter_nonMD_Morph(i)).ODI_AMPA_step_peak;
    ODI_N(i)=data(filter_nonMD_Morph(i)).ODI_NMDA_step_peak;

end
NMDA_silent = abs(ODI_A)<1 & abs(ODI_N)==1;

if ignore_missing_Nresp_cells
    ODI_A(NMDA_silent)=nan;
    ODI_N(NMDA_silent)=nan;
end
ODI_usil = nan(1,length(ODI_A));
ODI_usil(ODI_A>0) = (3-5.*ODI_N(ODI_A>0))./(3.*ODI_N(ODI_A>0)-5);
ODI_usil(ODI_A<0) = (5.*ODI_N(ODI_A<0)+3)./(3.*ODI_N(ODI_A<0)+5);

dom_AN_ratio(NMDA_silent)=nan; %% inf A/N cannot be processed
nondom_AN_ratio(NMDA_silent)=nan;

mon_cells = abs(ODI_A)==1;
mono_AN = dom_AN_ratio(mon_cells);
domi_AN = dom_AN_ratio(~mon_cells);
nondomi_AN = nondom_AN_ratio(~mon_cells);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'AMPAvsNMDA', 'Position', [1467,657,437,318]); 
fig_number = fig_number+1;
subplot(1,2,1)
plotSpread({1-abs(ODI_N),1-abs(ODI_A)},'categoryLabels',{'NMDA','AMPA'},'distributionColors',{'k','k'})
for  i = 1:length(filter_nonMD_Morph)
    plot([1,2],[1-abs(ODI_N(i)),1-abs(ODI_A(i))],'color',[0.5 0.5 0.5]); hold on
end
plot([1,2],nanmedian([1-abs(ODI_N(:)),1-abs(ODI_A(:))]),'k','LineWidth',3)
set(gca,'TickDir','out')
xticklabels({'NMDA','AMPA'})
ylabel('1-abs(ODI)')
ylim([-0.01 1.01]); 
[p,~,STATS]=signrank([(1-abs(ODI_A(:)))],[(1-abs(ODI_N(:)))],'tail', 'both');
title({'Wilcoxon signrank' ['(n=' num2str(sum(~isnan(ODI_A))) ', p=' num2str(p,2)] ...
    [ 'W=' num2str(STATS.signedrank) ')']})

subplot(1,2,2); cla
plotSpread({mono_AN,domi_AN,nondomi_AN},'categoryLabels',{'mon','dom','nondom'},'distributionColors',{'k','k','k'})
for  i = 1:length(domi_AN)
    if ~isnan(domi_AN(i))
        plot([2,3],[domi_AN(i),nondomi_AN(i)],'color',[0.5 0.5 0.5]); hold on
    end
end
plot([2,3],nanmedian([domi_AN;nondomi_AN],2),'k','LineWidth',3)
plot([0.7,1.3],nanmedian([mono_AN;mono_AN],2),'k','LineWidth',3)
plot(xlim, [1 1],'--k')
set(gca,'TickDir','out')
xticklabels({'mon','dom','non-dom'})
ylabel('AMPA/NMDA'); ylim([-0.1 inf])
[p_mon_dom,~,STATS1]=ranksum(mono_AN,domi_AN,'tail', 'both');
[p_dom_nondom,~,STATS2]=signrank(domi_AN,nondomi_AN,'tail', 'both');
title({['MWU (n=' num2str(sum(~isnan(mono_AN))) ' & ' num2str(sum(~isnan(domi_AN)))...
    ', p=' num2str(p_mon_dom,2) ','], [',U=' num2str(STATS1.ranksum) ')'],...
    'Wilcoxon signrank', ['(n= ' num2str(sum(~isnan(domi_AN))) ', p=' num2str(p_dom_nondom,2)...
    ', W=' num2str(STATS2.signedrank) ')']})

%% save figs
if savefigs
    mkdir([savefigs_location '\synaptic_mech_figs\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[savefigs_location '\synaptic_mech_figs\' fig_hand(i).Name],'svg')
        saveas(fig_hand(i),[savefigs_location '\synaptic_mech_figs\' fig_hand(i).Name],'fig')
    end
end