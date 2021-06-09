if savefigs
    close all
    clear fig_hand
    fig_number = 1;
end

%% Dependencies
% The following code relies on these dependencies which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/
%   BeeswarmPlot 
%   cbrewer 
%   suptitle
%   CircStat2012a
%   circ_corrclnp

%% structure_funciton_cell_morph_fig_cirstats
filter_Morph = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', 1, 'dLGN', 1, 'A_ramp_signal_QC', nan, ... 
    'Conf_QC', 1, 'TwoPTrans_QC', 1, 'CS_internal',nan,'Interneuron_morph',0,...
    'morph_QC',1));

data_ex_cell = data(22);

Fr_rotation_example = data_ex_cell.FlourRatio_data.MorphDensity3D_rotated_adjusted_Fratio(:,1);
Morph_rot_angles = data_ex_cell.FlourRatio_data.Morph_rot_angles;
Fr_radialMask_example = data_ex_cell.FlourRatio_data.Radial3D_adjusted_Fratio;

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'mFR_over_rotations_example', 'Position', [1,451,827,542]); 
fig_number = fig_number+1;
subplot(2,4,[1:4]);
plot(Morph_rot_angles,Fr_rotation_example,'k','LineWidth', 2); hold on
plot([1 360], [mean(Fr_rotation_example) mean(Fr_rotation_example)],'--k','LineWidth', 2); hold on
xlim([0 360]); ylim([-0.333 0.333])
xlabel('rotation (deg)'); 
ylabel('FR')
set(gca,'TickDir','out'); box off
title({'FR change when rotating morphology'; 'and mean over all rotations(--)'})

rot_cardinal = [find(Morph_rot_angles==0),find(Morph_rot_angles==90),find(Morph_rot_angles==180),find(Morph_rot_angles==270)];
for rot_i = 1:length(rot_cardinal)
    rotation_deg = rot_cardinal(rot_i);
    Fr_binned = data_ex_cell.FlourRatio_data.Binned_ratio_image3D_all(:,:,:,1);
    Morph_binned = full(data_ex_cell.FlourRatio_data.ThreeD_dendrite_density_rotated(:,:,:,rotation_deg));
    Morph_binned(Morph_binned==0) = nan;
    Fr_binned(isnan(Morph_binned)) = nan; % apply morph mask
    Fr_binned_zproj = squeeze(nanmean(Fr_binned,3));
    Morph_binned_zproj = squeeze(nanmean(Morph_binned,3));
    ax(rot_i+4) = subplot(2,4,rot_i+4); cla
    h = imagesc(Fr_binned_zproj); hold on; axis image
    scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
    set(h,'AlphaData',...
        (Morph_binned_zproj./max(Morph_binned_zproj(:)).*(1-0.5)+0.5)...
        .*~isnan(Morph_binned_zproj)...
        .*~isnan(Fr_binned_zproj));
%     set(gca,'Color',[0.85 0.85 0.85])
    caxis([-0.7 0.7])
    xlim([21-10.5 21+10.5]); ylim([21-10.5 21+10.5])
    xticks([]); yticks([]);
    colormap(ax(rot_i+4),flip(cbrewer('div','RdYlGn',100),1))
    pbaspect([1 1 1]); set(gca,'TickDir','out'); %colorbar
    title([num2str(Morph_rot_angles(rot_cardinal(rot_i))) '°'])
end

radii = data(filter_Morph(1)).FlourRatio_data.Sphere_radii(1:end);
Morph_rot_angles = data(filter_Morph(1)).FlourRatio_data.Morph_rot_angles(1:end);
% Morph_rot_angles(Morph_rot_angles>180) = Morph_rot_angles(Morph_rot_angles>180)-360;

% example cell background vs asymetry
clear angle_of_grad_axons angle_of_grad_dendrites axon_angle_arrow_180  dend_angle_arrow_180
Binned_ratio_image3D_150um_cell = data_ex_cell.FlourRatio_data.Binned_ratio_image3D_all(:,:,:,radii==150);
Binned_ratio_image3D_150um_cell = nanmean(Binned_ratio_image3D_150um_cell,3);
center_pix_xy = ceil(size(Binned_ratio_image3D_150um_cell,1)/2);

% background
[x,y] = ndgrid(1:size(Binned_ratio_image3D_150um_cell,1),1:size(Binned_ratio_image3D_150um_cell,1));
coef = regress(Binned_ratio_image3D_150um_cell(:),cat(2,x(:),y(:))-center_pix_xy); % center at [0,0]
angle_of_grad_axons = cart2pol(coef(2),coef(1)); % ranges from -pi to +pi (direction)

[axon_angle_arrow(1), axon_angle_arrow(2)] = pol2cart(angle_of_grad_axons,10);
axon_angle_arrow = axon_angle_arrow+center_pix_xy; % center on soma
[axon_angle_arrow_orth(1), axon_angle_arrow_orth(2)] = pol2cart(angle_of_grad_axons+pi,10);
axon_angle_arrow_orth = axon_angle_arrow_orth+center_pix_xy; % center on soma

% dendrites 
morph = data_ex_cell.FlourRatio_data.morphology_transformed_scalled([2,1,3],:)';
morph(:,1) = morph(:,1)-373/2; % center at 0,0
morph(:,2) = morph(:,2)-373/2;
scaling_factor = size(data_ex_cell.Fr_sphere_stacks.Stack_red)./size(data_ex_cell.FlourRatio_data.Binned_ratio_image3D_all(:,:,:,1));
morph = morph./scaling_factor;

dend_map = full(data_ex_cell.FlourRatio_data.TwoD_dendrite_density_rotated(:,:,1,radii==150));

elongation = data_ex_cell.morphology.elongation_ratio;
elongation_angle = data_ex_cell.morphology.elongation_angle;

[dend_elong_angle_arrow(1), dend_elong_angle_arrow(2)] = pol2cart(elongation_angle,10);
dend_elong_angle_arrow = dend_elong_angle_arrow+center_pix_xy;
[dend_elong_angle_arrow_orth(1), dend_elong_angle_arrow_orth(2)] = pol2cart(elongation_angle+pi,10);
dend_elong_angle_arrow_orth = dend_elong_angle_arrow_orth+center_pix_xy;

asym_ratio = data_ex_cell.morphology.asymetry_ratio;
centermass_angle = data_ex_cell.morphology.asymetry_angle;

[dend_asim_angle_arrow(1), dend_asim_angle_arrow(2)] = pol2cart(centermass_angle,10);
dend_asim_angle_arrow = dend_asim_angle_arrow+center_pix_xy;
[dend_asim_angle_arrow_orth(1), dend_asim_angle_arrow_orth(2)] = pol2cart(centermass_angle+pi,10);
dend_asim_angle_arrow_orth = dend_asim_angle_arrow_orth+center_pix_xy;

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'BGvsElongvAsym_example', 'Position', [1,38,420,423]); 
fig_number = fig_number+1;
ax(1) = subplot(1,2,1); cla
h = imagesc(Binned_ratio_image3D_150um_cell); hold on
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
colormap(ax(1),flip(cbrewer('div','RdYlGn',100),1))
set(h,'AlphaData',~isnan(Binned_ratio_image3D_150um_cell))
plot([axon_angle_arrow(1),axon_angle_arrow_orth(1)],[axon_angle_arrow(2),axon_angle_arrow_orth(2)],'b','LineWidth',2)
xticks([]); yticks([]); axis on image
caxis([-0.7 0.7]); box on
xlim([21-10.5 21+10.5]); ylim([21-10.5 21+10.5])
title(['angle: ' num2str(wrapTo180(rad2deg(wrapToPi(angle_of_grad_axons))),3)  '°'])

ax(2) = subplot(1,2,2); cla
h = imagesc(dend_map); hold on
scatter(data_ex_cell.FlourRatio_data.Binned_cell_pos(1),data_ex_cell.FlourRatio_data.Binned_cell_pos(2),[],'b','filled')
scatter(morph(:,1)+data_ex_cell.FlourRatio_data.Binned_cell_pos(1),...
    morph(:,2)+data_ex_cell.FlourRatio_data.Binned_cell_pos(2),1,'k','filled')
colormap(ax(2),flip(colormap(ax(2),'gray')))
plot([dend_elong_angle_arrow(1),dend_elong_angle_arrow_orth(1)],[dend_elong_angle_arrow(2),dend_elong_angle_arrow_orth(2)],'--b','LineWidth',2)
plot([dend_asim_angle_arrow(1),dend_asim_angle_arrow_orth(1)],[dend_asim_angle_arrow(2),dend_asim_angle_arrow_orth(2)],'b','LineWidth',2)

xticks([]); yticks([]); axis on image
set(ax(2),'Color',[1 1 1]); box on
xlim([21-10.5 21+10.5]); ylim([21-10.5 21+10.5])
title({['asym (b): ' num2str(wrapTo180(rad2deg(wrapToPi(centermass_angle+pi))),3) '°'];...
    ['elong (--b): ' num2str(wrapTo180(rad2deg(wrapToPi(elongation_angle))),3) '°']})

angle_dif_asym = circ_dist(angle_of_grad_axons,centermass_angle+pi); % difference in ori (not dir (thats what the *2 is for))
angle_dif_elong = circ_dist(angle_of_grad_axons,elongation_angle);
suptitle({['angle dif (asym): ' num2str(wrapTo180(rad2deg(wrapToPi(angle_dif_asym))),3)  '°'];...
    ['angle dif (elong): ' num2str(wrapTo180(rad2deg(wrapToPi(angle_dif_elong))),3)  '°']});

%% population analysis
clear angle_of_grad_axons angle_of_grad_dendrites axon_angle_arrow dend_angle_arrow
clear centermass_angle elongation_angle angle_dif_asym angle_dif_elong asym_ratio1 asym_ratio2 asym_ratio
clear centermass_angle1 centermass_angle2 mag_of_grad_axons Conf_slice_align_angle
showplots = 0;
for i = 1:length(filter_Morph)
    
    % dendrite morph
    morph = data(filter_Morph(i)).FlourRatio_data.morphology_transformed_scalled([2,1,3],:)';
    morph(:,1) = morph(:,1)-373/2; % center at 0,0
    morph(:,2) = morph(:,2)-373/2;
    scaling_factor = size(data(filter_Morph(i)).Fr_sphere_stacks.Stack_red)./size(data(filter_Morph(i)).FlourRatio_data.Binned_ratio_image2D_all);
    morph = morph./scaling_factor;
    
    dend_map = full(data(filter_Morph(i)).FlourRatio_data.TwoD_dendrite_density_rotated(:,:,1,radii==150));

    % background angle
    Binned_ratio_image3D_150um_cell = data(filter_Morph(i)).FlourRatio_data.Binned_ratio_image3D_all(:,:,:,radii==150);
    Binned_ratio_image3D_150um_cell = nanmean(Binned_ratio_image3D_150um_cell,3);
    center_pix_xy = ceil(size(Binned_ratio_image3D_150um_cell,1)/2);
    
    [x,y] = ndgrid(1:size(Binned_ratio_image3D_150um_cell,1),1:size(Binned_ratio_image3D_150um_cell,1));
    coef = regress(Binned_ratio_image3D_150um_cell(:),cat(2,x(:),y(:))-center_pix_xy); % center at [0,0]
    [background_grad_angle, background_grad_mag] = cart2pol(coef(2),coef(1));
    
    elongation(i) = data(filter_Morph(i)).morphology.elongation_ratio;
    elongation_angle(i) = wrapToPi(data(filter_Morph(i)).morphology.elongation_angle);
    asym_ratio(i) = data(filter_Morph(i)).morphology.asymetry_ratio;
    centermass_angle(i) = wrapToPi(data(filter_Morph(i)).morphology.asymetry_angle);  
    mag_of_grad_axons(i) = background_grad_mag;
    angle_of_grad_axons(i) = background_grad_angle;
    Conf_slice_align_angle(i) = wrapToPi(deg2rad(data(filter_Morph(i)).Conf_to_CCF_rotation));
    
end

% conf to ccf angle offset adjustment
elongation_angle_adj = elongation_angle+Conf_slice_align_angle;
centermass_angle_adj = centermass_angle+Conf_slice_align_angle;
angle_of_grad_axons_adj = angle_of_grad_axons+Conf_slice_align_angle;

% convert from dir to ori
elongation_angle = wrapToPi(elongation_angle.*2)./2;  % non adjusted for accf alignement
centermass_angle = wrapToPi(centermass_angle.*2)./2; 
angle_of_grad_axons = wrapToPi(angle_of_grad_axons.*2)./2; 

centermass_angle_adj_dir = wrapToPi(centermass_angle_adj); % dir not ori
angle_of_grad_axons_adj_dir = wrapToPi(angle_of_grad_axons_adj); 

elongation_angle_adj = wrapToPi(elongation_angle_adj.*2)./2; % adjusted and converted to ori
centermass_angle_adj = wrapToPi(centermass_angle_adj.*2)./2; 
angle_of_grad_axons_adj = wrapToPi(angle_of_grad_axons_adj.*2)./2; 


% calc angle diff
angle_dif_asym = circ_dist(angle_of_grad_axons,centermass_angle);
angle_dif_elong = circ_dist(angle_of_grad_axons,elongation_angle);

% convert from dir to ori (most values are already fine though)
angle_dif_asym = wrapToPi(angle_dif_asym.*2)./2; 
angle_dif_elong = wrapToPi(angle_dif_elong.*2)./2; 

figure;
scatter(elongation,asym_ratio,'b','filled')
xlim([0 1]); ylim([0 1])
xlabel('elong'); ylabel('asym')

%%
figure 
clear ax
subplot(1,3,1); polarhistogram(angle_of_grad_axons_adj_dir,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(1) = gca; title({'background FR '; 'gradiant dir'}); 
hold on; polarplot(repmat(circ_mean(angle_of_grad_axons_adj_dir'),1,2),(rlim),'--k','LineWidth',3)
subplot(1,3,2); polarhistogram(centermass_angle_adj_dir,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(2) = gca; title({'dendritic'; 'asym. dir'}); 
hold on; polarplot(repmat(circ_mean(centermass_angle_adj_dir'),1,2),(rlim),'--k','LineWidth',3)
subplot(1,3,3); polarhistogram(elongation_angle_adj,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(3) = gca; title({'dendritic'; 'elong. ori'}); 
hold on; polarplot(repmat(circ_mean(elongation_angle_adj'.*2)./2,1,2),(rlim),'--k','LineWidth',3)
for i = 1:2
    ax(i).ThetaLim=[-180 180];
    ax(i).ThetaTick=(-180:90:180);
    ax(i).ThetaTickLabels={'-90°','0°','90°','+/-180°'};
    ax(i).ThetaDir='clockwise';
    ax(i).RAxisLocation = 270;
end
ax(3).ThetaLim=[-90 90];
ax(3).ThetaTick=(-90:90:90);
ax(3).ThetaTickLabels={'-90°','0°','90°','+/-180°'};
ax(3).ThetaDir='clockwise';
ax(3).RAxisLocation = 270;


fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'BGvsElongvsAsym_distibution_ori', 'Position', [422,130,521,236]); 
fig_number = fig_number+1; 
clear ax
subplot(1,3,1); polarhistogram(angle_of_grad_axons_adj,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(1) = gca; 
title({'background FR '; 'gradiant ori'; ...
    ['mean: ' num2str(wrapTo180(rad2deg(circ_mean(angle_of_grad_axons_adj'.*2)./2)),3)]}); 
hold on; polarplot(repmat(circ_mean(angle_of_grad_axons_adj'.*2)./2,1,2),(rlim),'--k','LineWidth',3)

subplot(1,3,2); polarhistogram(centermass_angle_adj,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(2) = gca; 
title({'background FR '; 'asym. ori'; ...
    ['mean: ' num2str(wrapTo180(rad2deg(circ_mean(centermass_angle_adj'.*2)./2)),3)]}); 
hold on; polarplot(repmat(circ_mean(centermass_angle_adj'.*2)./2,1,2),(rlim),'--k','LineWidth',3)

subplot(1,3,3); polarhistogram(elongation_angle_adj,-pi:pi/8:pi,'FaceColor',[0.5 0.5 0.5]); ax(3) = gca; 
title({'background FR '; 'elong. ori'; ...
    ['mean: ' num2str(wrapTo180(rad2deg(circ_mean(elongation_angle_adj'.*2)./2)),3)]}); 
hold on; polarplot(repmat(circ_mean(elongation_angle_adj'.*2)./2,1,2),(rlim),'--k','LineWidth',3)
for i = 1:length(ax)
    ax(i).ThetaLim=[-90 90];
    ax(i).ThetaTick=(-90:90:90);
    ax(i).ThetaTickLabels={'-90°','0°','90°','+/-180°'};
    ax(i).ThetaDir='clockwise';
    ax(i).RAxisLocation = 270;
end

angle_of_grad_axons_centered = wrapToPi(angle_of_grad_axons_adj*2-circ_mean(angle_of_grad_axons_adj'*2))./2;
centermass_angle_centered = wrapToPi(centermass_angle_adj*2-circ_mean(centermass_angle_adj'*2))./2;
elongation_angle_centered = wrapToPi(elongation_angle_adj*2-circ_mean(elongation_angle_adj'*2))./2;

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'BG_Elong_Asym_corr', 'Position', [829,673,530,318]); 
fig_number = fig_number+1;
subplot(1,2,1)
scatter(angle_of_grad_axons_centered/pi,centermass_angle_centered/pi,[],'k','filled','MarkerFaceAlpha',0.8); axis equal square
xlabel({'bg FR ori.' '(mean subtracted deg)'});ylabel({'dend. asym. ori.' '(mean subtracted deg)'}); set(gca,'TickDir','out')
[r, corr_p]=circ_corrcc(angle_of_grad_axons_centered*2,centermass_angle_centered*2); % we care about corr because it looks how the two variables co-vary
title({['background vs asym']; ['n: ' num2str(length(angle_of_grad_axons_centered)) ', r: ' num2str(r,2) ', p: ' num2str(corr_p,2)]})
xticks([-0.5 0 +0.5]); xticklabels({'-90°' '0°' '+90°'})
yticks([-0.5 0 +0.5]); yticklabels({'-90°' '0°' '+90°'})

subplot(1,2,2)
scatter(angle_of_grad_axons_centered/pi,elongation_angle_centered/pi,[],'k','filled','MarkerFaceAlpha',0.8); axis equal square
xlabel({'bg FR ori.' '(mean subtracted deg)'});ylabel({'dend. elong. ori.' '(mean subtracted deg)'}); set(gca,'TickDir','out')
[r, corr_p]=circ_corrcc(angle_of_grad_axons_centered*2,elongation_angle_centered*2); % we care about corr because it looks how the two variables co-vary
title({['background vs elong']; ['n: ' num2str(length(angle_of_grad_axons_centered)) ', r: ' num2str(r,2) ', p: ' num2str(corr_p,2)]})
xticks([-0.5 0 +0.5]); xticklabels({'-90°' '0°' '+90°'})
yticks([-0.5 0 +0.5]); yticklabels({'-90°' '0°' '+90°'})


figure;
subplot(2,2,1)
scatter(angle_of_grad_axons_centered/pi,centermass_angle_centered/pi,[],'k','filled','MarkerFaceAlpha',0.8); axis equal square
xlabel('FR gradient (rad/pi)');ylabel('dend asym (rad/pi)'); set(gca,'TickDir','out')
[r, corr_p]=circ_corrcc(angle_of_grad_axons_centered*2,centermass_angle_centered*2); % we care about corr because it looks how the two variables co-vary
title({['background vs asym']; ['r: ' num2str(r,2) ', p: ' num2str(corr_p,2)]})

subplot(2,2,2)
scatter(angle_of_grad_axons_centered/pi,elongation_angle_centered/pi,[],'k','filled','MarkerFaceAlpha',0.8); axis equal square
xlabel('FR gradient (rad/pi)');ylabel('dend elong (rad/pi)');set(gca,'TickDir','out')
[r, corr_p]=circ_corrcc(angle_of_grad_axons_centered*2,elongation_angle_centered*2); % we care about corr because it looks how the two variables co-vary
title({['background vs elong']; ['r: ' num2str(r,2) ', p: ' num2str(corr_p,2)]})

subplot(2,2,3)
AP_pos = cellfun(@(x) x(3), {data(filter_Morph(:)).ccfv3_pos});
scatter(AP_pos,angle_of_grad_axons_centered/pi,[],'k','filled','MarkerFaceAlpha',0.8); axis equal square
xlabel('AP pos');ylabel('background grad (rad/pi)');set(gca,'TickDir','out')
[r, corr_p]=circ_corrcl(angle_of_grad_axons_centered*2,AP_pos); % we care about corr because it looks how the two variables co-vary
title({['AP pos vs background grad']; ['r: ' num2str(r,2) ', p: ' num2str(corr_p,2)]})

subplot(2,2,4)
scatter(AP_pos,elongation_angle_centered/pi,[],'k','filled','MarkerFaceAlpha',0.8); axis equal square
xlabel('AP pos');ylabel('dend elong (rad/pi)');set(gca,'TickDir','out')
[r, corr_p]=circ_corrcl(elongation_angle_centered*2,AP_pos); % we care about corr because it looks how the two variables co-vary
title({['AP pos vs elongation']; ['r: ' num2str(r,2) ', p: ' num2str(corr_p,2)]})

%%
filter_Morph_temp = find(QC_cell_filter(data, your_data_dir, ...
    'transduction_QC', 1, 'dLGN', 1, 'A_ramp_signal_QC', nan, ... 
    'Conf_QC', 1, 'TwoPTrans_QC', 1, 'CS_internal',nan,'Interneuron_morph',0,...
    'morph_QC',1));
filter_Morph = filter_Morph_temp;

perc_morph_cut_norot = cell2mat(cellfun(@(x) x.Binned_ratio_images3D_percent_dendrites_excluded(1,1), {data((filter_Morph)).FlourRatio_data},'UniformOutput',false));
perc_morph_cut_allrot_min = cell2mat(cellfun(@(x) min(x.Binned_ratio_images3D_percent_dendrites_excluded(:,1)), {data((filter_Morph)).FlourRatio_data},'UniformOutput',false));
perc_morph_cut_allrot_max = cell2mat(cellfun(@(x) max(x.Binned_ratio_images3D_percent_dendrites_excluded(:,1)), {data((filter_Morph)).FlourRatio_data},'UniformOutput',false));
perc_morph_cut_allrot_dif = perc_morph_cut_allrot_max-perc_morph_cut_allrot_min;
perc_morph_cut_norot_cut=30;
perc_morph_cut_allrot_dif_cut=10;
filter_Morph([perc_morph_cut_norot>perc_morph_cut_norot_cut] | [perc_morph_cut_allrot_dif>perc_morph_cut_allrot_dif_cut]) = []; % remove cells with too many dendrites cutt off

figure; set(gcf,'color','w')
scatter(perc_morph_cut_norot([perc_morph_cut_norot>perc_morph_cut_norot_cut] | [perc_morph_cut_allrot_dif>perc_morph_cut_allrot_dif_cut]),...
    perc_morph_cut_allrot_dif([perc_morph_cut_norot>perc_morph_cut_norot_cut] | [perc_morph_cut_allrot_dif>perc_morph_cut_allrot_dif_cut]),[],'b','filled',...
    'MarkerFaceAlpha',0.4); hold on
scatter(perc_morph_cut_norot(~[perc_morph_cut_norot>perc_morph_cut_norot_cut] & ~[perc_morph_cut_allrot_dif>perc_morph_cut_allrot_dif_cut]),...
    perc_morph_cut_allrot_dif(~[perc_morph_cut_norot>perc_morph_cut_norot_cut] & ~[perc_morph_cut_allrot_dif>perc_morph_cut_allrot_dif_cut]),[],'b','filled');
plot([min(xlim) max(xlim)],[10 10],'--k')
plot([30 30],[min(ylim) max(ylim)],'--k')
set(gca,'TickDir','out')
xlabel('% dend. cut at oritinal orientation')
ylabel('max(delta(%)) dend. cut with morph rotation')
title({['thresholds for exlcuding morphologies'] ['based on dendrites cut off']})

% now looking at the bias of axon sampling by the dendrites using
% correlation of even sampling vs dendritic sampling at different rotations
% rotate morph
temp = cellfun(@(x) x.MorphDensity3D_rotated_adjusted_Fratio(:,1), {data((filter_Morph)).FlourRatio_data},'UniformOutput',false);
ThreeD_MorphDens_rot = cat(2,temp{:})';
temp = cellfun(@(x) x.Radial3D_adjusted_Fratio(:,1), {data((filter_Morph)).FlourRatio_data},'UniformOutput',false);
Radial_MorphDens_rot = cat(2,temp{:})';
temp = cellfun(@(x) x.MeanRadial3D_adjusted_Fratio(:,1), {data((filter_Morph)).FlourRatio_data},'UniformOutput',false);
MeanRadial_MorphDens_rot = cat(2,temp{:})';
Morph_rot_angles = data_ex_cell.FlourRatio_data.Morph_rot_angles;

mFR_mrFR_dif=tanh(atanh(ThreeD_MorphDens_rot)-atanh(MeanRadial_MorphDens_rot));
mFR_rFR_dif=tanh(atanh(ThreeD_MorphDens_rot)-atanh(Radial_MorphDens_rot));
mFR_mFR_dif=tanh(atanh(ThreeD_MorphDens_rot)-atanh(mean(ThreeD_MorphDens_rot,2))); % difference to mean of all rotations
mFR_mFR_dif_max = max(abs(mFR_mFR_dif),[],2);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'rotation_vs_monoc', 'Position', [1354,567,560,420]); 
fig_number = fig_number+1;
temp1 = ThreeD_MorphDens_rot;
subplot(2,2,1)
histogram(temp1(:,1), [-1:2/8:1],'FaceColor',[0.5 0.5 0.5])
xlabel('mFR at 0° rotation'); ylabel('cell count'); xlim([-1 1])
box off
pbaspect([2 1 1])
title(['mFR distribution (n: ' num2str(size(temp1,1)) ')'])
subplot(2,2,3)
histogram(abs(temp1(:,1)), [0:2/16:1],'FaceColor',[0.5 0.5 0.5])
xlabel('|mFR| at 0° rotation'); ylabel('cell count'); xlim([0 1])
box off; hold on
plot([median(abs(temp1(:,1))),median(abs(temp1(:,1)))],ylim,'--k','LineWidth',3)
pbaspect([2 1 1])
subplot(2,2,[2,4])
temp2 = repmat(deg2rad(Morph_rot_angles),size(temp1,1),1);
for i = 1:size(temp1,1)
    polarplot(temp2(i,:),abs(temp1(i,:)),'color',[0.5 0.5 0.5]); hold all
end
polarplot(temp2(1,:),median(abs(temp1),1),'--k','LineWidth',3)
thetaticks(0:90:360)
rlim([0 1])
thetaticklabels({'0°','90°','+/-180°','-90°'})
a = gca;
a.FontSize = 10;
a.RAxis.FontSize = 10;
a.ThetaDir='clockwise';
a.RAxisLocation = 270;
[rho, pval] = circ_corrclnp(wrapToPi(temp2(:)*2),abs(temp1(:))); 
title({'|mFR| over morph rotations';...
    ['np cl corr R: ' num2str(rho,2) ', p: ' num2str(pval,2)]})


percentile_grad = 75;
fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'max_FR_change_hist', 'Position', [1009,94,289,275]); 
fig_number = fig_number+1;
histogram(mFR_mFR_dif_max,10,'FaceColor',[0.5 0.5 0.5]); hold on
xlim([0 inf]); xlabel('magnitude of background gradient'); ylabel('cell count')
plot([prctile(mFR_mFR_dif_max,percentile_grad), prctile(mFR_mFR_dif_max,percentile_grad)], ylim,'--k','LineWidth',3)
box off
axis square
title(['threshold at ' num2str(percentile_grad) '%: ' num2str(prctile(mFR_mFR_dif_max,percentile_grad),2)])

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'rotation_vs_monoc_subgroup', 'Position', [1354,58,560,420]); 
fig_number = fig_number+1;
temp1 = ThreeD_MorphDens_rot(mFR_mFR_dif_max>(prctile(mFR_mFR_dif_max,percentile_grad)),:);
subplot(2,2,1)
histogram(temp1(:,1), [-1:2/8:1],'FaceColor',[0.5 0.5 0.5])
xlabel('mFR at 0° rotation'); ylabel('cell count'); xlim([-1 1])
set(gca,'TickDir','out'); box off
title(['mFR distribution (n: ' num2str(size(temp1,1)) ')'])
pbaspect([2 1 1])
subplot(2,2,3)
histogram(abs(temp1(:,1)), [0:2/16:1],'FaceColor',[0.5 0.5 0.5])
xlabel('|mFR| at 0° rotation'); ylabel('cell count'); xlim([0 1])
set(gca,'TickDir','out'); box off; hold on
plot([median(abs(temp1(:,1))),median(abs(temp1(:,1)))],ylim,'--k','LineWidth',3)
pbaspect([2 1 1])
subplot(2,2,[2,4])
temp2 = repmat(deg2rad(Morph_rot_angles),size(temp1,1),1);
for i = 1:size(temp1,1)
    polarplot(temp2(i,:),abs(temp1(i,:)),'color',[0.5 0.5 0.5]); hold all
end
polarplot(temp2(1,:),median(abs(temp1),1),'--k','LineWidth',3)
thetaticks(0:90:360)
rlim([0 1])
thetaticklabels({'0°','90°','+/-180°','-90°'})
a = gca;
a.FontSize = 10;
a.RAxis.FontSize = 10;
a.ThetaDir='clockwise';
a.RAxisLocation = 270;
[rho, pval] = circ_corrclnp(wrapToPi(temp2(:)*2),abs(temp1(:))); 
title({'|mFR| over morph rotations'; 'only including a subset of cells';...
    ['np cl corr R: ' num2str(rho,2) ', p: ' num2str(pval,2)]}) 


%% save figs
if savefigs
    mkdir([savefigs_location '\monoc_via_morph\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[savefigs_location '\monoc_via_morph\' fig_hand(i).Name],'svg')
        saveas(fig_hand(i),[savefigs_location '\monoc_via_morph\' fig_hand(i).Name],'fig')
    end
end