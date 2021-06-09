if savefigs
    close all
    clear fig_hand
    fig_number = 1;
end

%% Filter for paired AMPA/NMDA, also filter for eye injections
am_nm = QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', 1,...
    'dLGN', 1, ...
    'A_ramp_signal_QC', 1,...
    'N_ramp_signal_QC', 1, ...
    'AtoN_Rschange_QC',1,...
    'CS_internal', 1);

am_nm1=QC_cell_filter(data, your_data_dir, ...
    'brain_contra_ipsi',1,...
    'transduction_QC', 1, ...
    'dLGN', 1, ...
    'A_ramp_signal_QC', 1, ...
    'N_ramp_signal_QC', 1, ...
    'AtoN_Rschange_QC',1,...
    'CS_internal', 1);

am_nm2=QC_cell_filter(data, your_data_dir,...
    'brain_contra_ipsi',0,...
    'transduction_QC', 1,...
    'dLGN', 1, ...
    'A_ramp_signal_QC', 1,...
    'N_ramp_signal_QC', 1,...
    'AtoN_Rschange_QC',1,...
    'CS_internal', 1);

%% Read out odi and ctaegorized based on ODI into contra, ipsi, bino, contra silent and ipsi silent to show example traces
%ODI AMPA and NMDA
odi_am_nm=[vertcat(data(am_nm).ODI_AMPA_step_peak) vertcat(data(am_nm).ODI_NMDA_step_peak)];
%conta
cat1= odi_am_nm(:,1)==1 & odi_am_nm(:,2)==1;
%ipsi
cat2= odi_am_nm(:,1)==-1 & odi_am_nm(:,2)==-1;
%bino
cat3= odi_am_nm(:,1)~=1  & odi_am_nm(:,2)~=1  & odi_am_nm(:,1)~=-1  & odi_am_nm(:,2)~=-1;
%ipsi silent
cat4= odi_am_nm(:,1)==1  & odi_am_nm(:,2)<1  & odi_am_nm(:,2)>0;
%contra silent
cat5= odi_am_nm(:,1)==-1  & odi_am_nm(:,2)>-1  & odi_am_nm(:,2)<0;

all_cat=cat1+cat2+cat3+cat4+cat5;

%% ODI distributions and category plot
colorscale_ODI = (cbrewer('div','RdYlBu',13)); %% for 6 bins

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'ODI_category_histogram', 'Position', [200, 600, 1000, 300]); 
fig_number = fig_number+1;

%Step AMPA/NMDA
subplot(1,3,1);
b=bar([1 2 3 4 5],[sum(cat1) sum(cat2) sum(cat3)...
    sum(cat4) sum(cat5)]);
box off;ylabel('Count');
b.FaceColor = 'flat';ylim([0 60]);yticks([0:30:60])
b.CData(1,:) = [0.7 0.7 0.7];b.CData(2,:) = [0.7 0.7 0.7]; b.CData(3,:) = [1 1 1]; b.CData(4,:) = [1 1 1];b.CData(5,:) = [1 1 1];
yticks([0:30:60]); xticklabels({'contra','ipsi','bino','ipsi silent','contra silent'});set(gca,'FontSize',10)
xtickangle(45);set(gca,'FontSize',10);

clear ax
ax(1) = subplot(1,3,2);
binos = abs(odi_am_nm(:,1))<1; 
binobincounts = histcounts(odi_am_nm(find(binos),1),[-1:2/6:1]); 
binscounts = [sum(odi_am_nm(:,1)==-1), 0, binobincounts, 0, sum(odi_am_nm(:,1)==1)];
b = bar(binscounts,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title(['AMPA, n=' num2str(length(odi_am_nm(:,1)))])

ax(2) = subplot(1,3,3);
binos = abs(odi_am_nm(:,2))<1; 
binobincounts = histcounts(odi_am_nm(find(binos),2),[-1:2/6:1]); 
binscounts = [sum(odi_am_nm(:,2)==-1), 0, binobincounts, 0, sum(odi_am_nm(:,2)==1)];
b = bar(binscounts,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title(['NMDA, n=' num2str(length(odi_am_nm(:,2)))])

for i = 1:length(ax)
    set(ax(i),'XTick', [1,2.5,5.5,8.5 10])
    set(ax(i),'XTickLabels', {'-1','>-1','0','<1','1'})
    set(ax(i),'box','off')
    ax(i).YLabel.String='count';
    set(ax(i),'TickDir','out')
    ax(i).XLabel.String='ODI';
end

%% save figs
if savefigs
    mkdir([savefigs_location '\Categorical_vs_continouse_ODI\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[savefigs_location '\Categorical_vs_continouse_ODI\' fig_hand(i).Name],'svg')
        saveas(fig_hand(i),[savefigs_location '\Categorical_vs_continouse_ODI\' fig_hand(i).Name],'fig')
    end
end