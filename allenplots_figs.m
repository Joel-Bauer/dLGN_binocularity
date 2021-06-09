if savefigs
    close all
    clear fig_hand
    fig_number = 1;
end

%% Dependencies
% The following code relies on these dependencies which can be downloaded from
% https://de.mathworks.com/matlabcentral/fileexchange/
%   distinguishable_colors 
%   BeeswarmPlot 
%   cbrewer 

%% select cells for ODI alignment with retinotopy 
filter_ACCF_retino = find(QC_cell_filter(data, your_data_dir,...
    'transduction_QC', 1, ...
    'dLGN', 1, ...
    'A_ramp_signal_QC', 1, ...
    'N_ramp_signal_QC', nan,...
    'Conf_QC', nan, ...
    'CCF_QC', 1, ...
    'CS_internal', 1,...
    'Interneuron_morph',0));

ODI_temp = {data(filter_ACCF_retino).ODI_AMPA_step_peak};
ODI_temp = [ODI_temp{:}];

ccfpos_temp = {data(filter_ACCF_retino).ccfv3_pos};
ccfpos_temp = reshape([ccfpos_temp{:}],3,[])';

ccfpos_temp2 = cat(2,ccfpos_temp,ODI_temp'); % ccfpos is in RSA
ccfpos_temp2(any(arrayfun(@(x) isnan(x),ccfpos_temp2),2),:) = [];

shape_sampled_ODI=alphaShape(ccfpos_temp2(:,3),ccfpos_temp2(:,1),ccfpos_temp2(:,2),inf); % alphashape was 30

colorscale_ODI = (cbrewer('div','RdYlBu',13)); %% for 6 bins

figure; set(gcf,'color','w')
plot(shape_dLGN,'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.1,'EdgeColor','none'); hold on
plot(shape_sampled_ODI,'FaceColor',[0 0 1],'FaceAlpha',0.1,'EdgeColor','none')
plot(shape_sampled_ret,'FaceColor',[1 0 0],'FaceAlpha',0.1,'EdgeColor','none')
plot(shape_sampled_overlap,'FaceColor',[0 1 0],'FaceAlpha',0.5,'EdgeColor','none')
axis equal
xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');  
ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right'); 
zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
caxis([-20 40]);
axis on vis3d 
box on
grid off
view([0 -1 0]);
axis manual
title({'Sampled area';' (b: ODI, r: Retinotopy, g: Overlap)'})

%% Dräger defined binocular visuotopic area
% crude estimation based on Dräger 1979.
binocular_boarder = [0,	-40;...
    0,	-20;...
    20,	0;...
    40,	20;...
    60,	40;...
    80,	50;...
    100, 60];

%% interp parameters
interp_gaussian_std = 7.5;
interp_sphereR_threshold = 5;

%% interpolation of eye preference (in the code this is called ODI because i couldnt be bothered to change it)
x_ODI = ccfpos_temp2(:,3);

xyz = ccfpos_temp2(:,[3,1,2]);
v=ccfpos_temp2(:,4);
v=(v>0)*2-1; % eye preference

[xq,yq,zq]=meshgrid([min(shape_dLGN.Points(:,1)):1:max(shape_dLGN.Points(:,1))],...
    [750:1:max(shape_dLGN.Points(:,2))],...
    [min(shape_dLGN.Points(:,3)):1:max(shape_dLGN.Points(:,3))]); % shape is ARS 

% full interpolation + extrapolation
vq_inLGN_odi = nan(size(xq));
for i = 1:size(xq,1) %a
    for ii = 1:size(xq,2) %r,
        for iii = 1:size(xq,3) %s
            if sum(vecnorm(xyz - [xq(i,ii,iii) yq(i,ii,iii) zq(i,ii,iii)], 2, 2)<(interp_sphereR_threshold*interp_gaussian_std))>0
                if inShape(shape_dLGN,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    idx_points = find(vecnorm(xyz - [xq(i,ii,iii) yq(i,ii,iii) zq(i,ii,iii)], 2, 2)<(interp_sphereR_threshold*interp_gaussian_std));
                    point_std_dist = vecnorm(xyz(idx_points,:) - [xq(i,ii,iii) yq(i,ii,iii) zq(i,ii,iii)], 2, 2)/interp_gaussian_std;
                    weightings = exp(-(point_std_dist.^2)/2);
                    weightings = weightings./max(weightings);
                    vq_inLGN_odi(i,ii,iii)=sum(v(idx_points).*weightings)/sum(weightings);
                end
            end
        end
    end
end

% interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_space_odi = vq_inLGN_odi;
vq_APsampled_space_odi(xq<min(x_ODI)) = nan;
vq_APsampled_space_odi(xq>max(x_ODI)) = nan;

% only interpolation
vq_sampled_space_odi = vq_inLGN_odi;
for i = 1:size(vq_inLGN_odi,1) %a
    for ii = 1:size(vq_inLGN_odi,2) %r,
        for iii = 1:size(vq_inLGN_odi,3) %s
            if ~isnan(vq_sampled_space_odi(i,ii,iii))
                if ~inShape(shape_sampled_ODI,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_sampled_space_odi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% partial interpolation within overlaping ODI & retino sampled area
vq_overlapping_sampled_space_odi = vq_inLGN_odi;
for i = 1:size(vq_inLGN_odi,1) %a
    for ii = 1:size(vq_inLGN_odi,2) %r,
        for iii = 1:size(vq_inLGN_odi,3) %s
            if ~isnan(vq_overlapping_sampled_space_odi(i,ii,iii))
                if ~inShape(shape_sampled_overlap,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_overlapping_sampled_space_odi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% dLGN SHELL interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_shell_odi = vq_inLGN_odi;
vq_APsampled_shell_odi(xq<min(x_ODI)) = nan;
vq_APsampled_shell_odi(xq>max(x_ODI)) = nan;
for i = 1:size(vq_inLGN_odi,1) %a
    for ii = 1:size(vq_inLGN_odi,2) %r,
        for iii = 1:size(vq_inLGN_odi,3) %s
            if ~isnan(vq_APsampled_shell_odi(i,ii,iii))
                if ~inShape(shape_shell ,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_APsampled_shell_odi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% dLGN CORE contra interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_coreC_odi = vq_inLGN_odi;
vq_APsampled_coreC_odi(xq<min(x_ODI)) = nan;
vq_APsampled_coreC_odi(xq>max(x_ODI)) = nan;
for i = 1:size(vq_inLGN_odi,1) %a
    for ii = 1:size(vq_inLGN_odi,2) %r,
        for iii = 1:size(vq_inLGN_odi,3) %s
            if ~isnan(vq_APsampled_coreC_odi(i,ii,iii))
                if ~inShape(shape_contra_core,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_APsampled_coreC_odi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% dLGN CORE ipsi interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_coreI_odi = vq_inLGN_odi;
vq_APsampled_coreI_odi(xq<min(x_ODI)) = nan;
vq_APsampled_coreI_odi(xq>max(x_ODI)) = nan;
for i = 1:size(vq_inLGN_odi,1) %a
    for ii = 1:size(vq_inLGN_odi,2) %r,
        for iii = 1:size(vq_inLGN_odi,3) %s
            if ~isnan(vq_APsampled_coreI_odi(i,ii,iii))
                if ~inShape(shape_ipsi,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_APsampled_coreI_odi(i,ii,iii)=nan;
                end
            end
        end
    end
end

%% interpolation of Binocularity 1-|ODI| 
xyz = ccfpos_temp2(:,[3,1,2]);
v=ccfpos_temp2(:,4);
v=1-abs(v);

[xq,yq,zq]=meshgrid([min(shape_dLGN.Points(:,1)):1:max(shape_dLGN.Points(:,1))],...
    [750:1:max(shape_dLGN.Points(:,2))],...
    [min(shape_dLGN.Points(:,3)):1:max(shape_dLGN.Points(:,3))]); % shape is ARS 

% full interpolation + extrapolation
vq_inLGN_absodi = nan(size(xq));
for i = 1:size(xq,1) %a
    for ii = 1:size(xq,2) %r,
        for iii = 1:size(xq,3) %s
            if sum(vecnorm(xyz - [xq(i,ii,iii) yq(i,ii,iii) zq(i,ii,iii)], 2, 2)<(interp_sphereR_threshold*interp_gaussian_std))>0
                if inShape(shape_dLGN,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    idx_points = find(vecnorm(xyz - [xq(i,ii,iii) yq(i,ii,iii) zq(i,ii,iii)], 2, 2)<(interp_sphereR_threshold*interp_gaussian_std));
                    point_std_dist = vecnorm(xyz(idx_points,:) - [xq(i,ii,iii) yq(i,ii,iii) zq(i,ii,iii)], 2, 2)/interp_gaussian_std;
                    weightings = exp(-(point_std_dist.^2)/2);
                    weightings = weightings./max(weightings);
                    vq_inLGN_absodi(i,ii,iii)=sum(v(idx_points).*weightings)/sum(weightings);
                end
            end
        end
    end
end

% interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_space_absodi = vq_inLGN_absodi;
vq_APsampled_space_absodi(xq<min(x_ODI)) = nan;
vq_APsampled_space_absodi(xq>max(x_ODI)) = nan;

% only interpolation
vq_sampled_space_absodi = vq_inLGN_absodi;
for i = 1:size(vq_inLGN_absodi,1) %a
    for ii = 1:size(vq_inLGN_absodi,2) %r,
        for iii = 1:size(vq_inLGN_absodi,3) %s
            if ~isnan(vq_sampled_space_absodi(i,ii,iii))
                if ~inShape(shape_sampled_ODI,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_sampled_space_absodi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% partial interpolation within overlaping ODI & retino sampled area
vq_overlapping_sampled_space_absodi = vq_inLGN_absodi;
for i = 1:size(vq_inLGN_absodi,1) %a
    for ii = 1:size(vq_inLGN_absodi,2) %r,
        for iii = 1:size(vq_inLGN_absodi,3) %s
            if ~isnan(vq_overlapping_sampled_space_absodi(i,ii,iii))
                if ~inShape(shape_sampled_overlap,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_overlapping_sampled_space_absodi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% dLGN SHELL interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_shell_absodi = vq_inLGN_absodi;
vq_APsampled_shell_absodi(xq<min(x_ODI)) = nan;
vq_APsampled_shell_absodi(xq>max(x_ODI)) = nan;
for i = 1:size(vq_inLGN_absodi,1) %a
    for ii = 1:size(vq_inLGN_absodi,2) %r,
        for iii = 1:size(vq_inLGN_absodi,3) %s
            if ~isnan(vq_APsampled_shell_absodi(i,ii,iii))
                if ~inShape(shape_shell ,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_APsampled_shell_absodi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% dLGN CORE contra interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_coreC_absodi = vq_inLGN_absodi;
vq_APsampled_coreC_absodi(xq<min(x_ODI)) = nan;
vq_APsampled_coreC_absodi(xq>max(x_ODI)) = nan;
for i = 1:size(vq_inLGN_absodi,1) %a
    for ii = 1:size(vq_inLGN_absodi,2) %r,
        for iii = 1:size(vq_inLGN_absodi,3) %s
            if ~isnan(vq_APsampled_coreC_absodi(i,ii,iii))
                if ~inShape(shape_contra_core,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_APsampled_coreC_absodi(i,ii,iii)=nan;
                end
            end
        end
    end
end

% dLGN CORE ipsi interpolation + partial extrapolation within sampled AP bounds
vq_APsampled_coreI_absodi = vq_inLGN_absodi;
vq_APsampled_coreI_absodi(xq<min(x_ODI)) = nan;
vq_APsampled_coreI_absodi(xq>max(x_ODI)) = nan;
for i = 1:size(vq_inLGN_absodi,1) %a
    for ii = 1:size(vq_inLGN_absodi,2) %r,
        for iii = 1:size(vq_inLGN_absodi,3) %s
            if ~isnan(vq_APsampled_coreI_absodi(i,ii,iii))
                if ~inShape(shape_ipsi,xq(i,ii,iii),yq(i,ii,iii),zq(i,ii,iii)) % SAR to ars: 2,3,1 or y, z, x
                    vq_APsampled_coreI_absodi(i,ii,iii)=nan;
                end
            end
        end
    end
end

%% Indictation of Dräger deffined binocular visuotopic area in dLGN

idx_retmap_inbinoarea=inpolygon(vq_APsampled_space_azimuth,vq_APsampled_space_elevation,...
    [binocular_boarder(:,1);flip(binocular_boarder(:,1)*-1)],...
    [binocular_boarder(:,2);flip(binocular_boarder(:,2))]);

idx_retmap_inbinoarea = single(idx_retmap_inbinoarea);
idx_retmap_inbinoarea(isnan(vq_APsampled_space_azimuth))=nan;

figure;
set(gcf,'Renderer','opengl')
ax(1)=subplot(1,1,1); hold on
h = slice(xq,...
    yq,...
    zq,...
    idx_retmap_inbinoarea,... % vq_inLGN,... % 
    [min(shape_dLGN.Points(:,1)):10:max(shape_dLGN.Points(:,1))],...
    [750:10:max(shape_dLGN.Points(:,2))],...
    [min(shape_dLGN.Points(:,3)):10:max(shape_dLGN.Points(:,3))]); % is plotted in terms of ARS for ease of navigation so RSA to ARS: z,x,
alpha('color')
set(h,'EdgeColor','none','FaceAlpha',0.3);
set(gcf,'color','w')
plot(shape_dLGN,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeColor','none')
axis equal
xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');  
ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right'); 
zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
caxis([-1 1]); colormap(cbrewer('div','RdYlBu',100))
axis on vis3d manual
box on
grid off
view([30 35]);
title('ODI interpolation')
ax(1).Position(2:3) = ax(1).Position(2:3)*2/3;

%% ODI overview
fig_hand(fig_number) = figure; 
set(fig_hand(fig_number), 'Name', 'ODI_3D_cells', 'Position', [1,31,1920,973]);
fig_number = fig_number+1;
clear ax
ax(1)=subplot(1,2,1); hold on
plot(shape_contra,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeColor','none')% plot LGN shape
scatter3(ccfpos_temp(:,3),ccfpos_temp(:,1),ccfpos_temp(:,2),[],ODI_temp,'filled'); colormap(cbrewer('div','RdYlBu',100))
axis equal
xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');  
ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right'); 
zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
axis on vis3d manual
box on
view([30 35]);
title(['ODI, n = ' num2str(length(ODI_temp)) ' cells'])
ax(1).Position(2:3) = ax(1).Position(2:3)*2/3;

ax(2)=subplot(1,2,2); hold on
h = slice(xq,...
    yq,...
    zq,...
    vq_APsampled_space_odi,... % vq_inLGN,... % 
    [min(shape_dLGN.Points(:,1)):10:max(shape_dLGN.Points(:,1))],...
    [750:10:max(shape_dLGN.Points(:,2))],...
    [min(shape_dLGN.Points(:,3)):10:max(shape_dLGN.Points(:,3))]); % is plotted in terms of ARS for ease of navigation so RSA to ARS: z,x,
alpha('color')
set(h,'EdgeColor','none','FaceAlpha',0.3);
set(gcf,'color','w')
plot(shape_dLGN,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.2,'EdgeColor','none')
axis equal
xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');  
ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right'); 
zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
caxis([-1 1]); colormap(cbrewer('div','RdYlBu',100))
axis on vis3d manual
box on
grid off
view([30 35]);
title('ODI interpolation')
ax(2).Position(2:3) = ax(2).Position(2:3)*2/3;

Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(gcf, 'StoreTheLink', Link);

%% coronal section overview

slice_number = 7;
slice_stepsize = floor((max(shape_dLGN.Points(:,1))-min(shape_dLGN.Points(:,1)))/slice_number);
apslices = [min(shape_dLGN.Points(:,1)):slice_stepsize:max(shape_dLGN.Points(:,1))];


x_ret = shape_sampled_ret.Points(:,1);
x_ODI = shape_sampled_ODI.Points(:,1);

apslices(apslices<min([x_ODI;x_ret])|apslices>max([x_ODI;x_ret]))=[];
% apslices = unique(x_ret);

clear slice_shape_points slice_shape_points_ipsi dLGN_boundary dLGN_boundary_ipsi slice_shape_points_shell dLGN_boundary_shell
for slice_n = 1:length(apslices)
    % make dLGN outline
    slice_shape_points{slice_n} = shape_dLGN.Points([shape_dLGN.Points(:,1)>(apslices(slice_n)-2) ...
        & shape_dLGN.Points(:,1)<(apslices(slice_n)+2) ...
        & shape_dLGN.Points(:,2)>500],:);
    dLGN_boundary{slice_n}=boundary(slice_shape_points{slice_n}(:,2),slice_shape_points{slice_n}(:,3),1);
    % make dLGN ipsi outline
    slice_shape_points_ipsi{slice_n} = shape_ipsi.Points([shape_ipsi.Points(:,1)>(apslices(slice_n)-2) ...
        & shape_ipsi.Points(:,1)<(apslices(slice_n)+2) ...
        & shape_ipsi.Points(:,2)>500],:);
    dLGN_boundary_ipsi{slice_n}=boundary(slice_shape_points_ipsi{slice_n}(:,2),slice_shape_points_ipsi{slice_n}(:,3),1);
    % make dLGN shell outline
    slice_shape_points_shell{slice_n} = shape_shell.Points([shape_shell.Points(:,1)>(apslices(slice_n)-2) ...
        & shape_shell.Points(:,1)<(apslices(slice_n)+2) ...
        & shape_shell.Points(:,2)>500],:);
    dLGN_boundary_shell{slice_n}=boundary(slice_shape_points_shell{slice_n}(:,2),slice_shape_points_shell{slice_n}(:,3),1);



end

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'ODI_ret_slices', 'Position', [601,50,828,921]);
fig_number = fig_number+1;
for slice_n = 1:length(apslices)
    displayed_slice = slice_n;
    clear ax
    ax(1)=subplot(length(apslices),4,1+(4*(slice_n-1))); cla
    hold on
    h = slice(xq,...
        yq,...
        zq,...
        vq_APsampled_space_odi,... % vq_inLGN_odi,...% 
        [apslices(slice_n)],...
        [],...
        []); % is plotted in terms of ARS for ease of navigation so RSA to ARS: z,x,
    alpha('color')
    set(h,'EdgeColor','none','FaceAlpha',1);
    set(gcf,'color','w')
        % plot dLGN outline
        plot3(ones(size(dLGN_boundary{slice_n})).*apslices(slice_n),slice_shape_points{slice_n}(dLGN_boundary{slice_n},2),...
            slice_shape_points{slice_n}(dLGN_boundary{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_ipsi{slice_n})).*apslices(slice_n),slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},2),...
            slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_shell{slice_n})).*apslices(slice_n),slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},2),...
            slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},3),'k','LineWidth',2)
    axis equal
    xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');
    ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right');
    zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
    caxis([-1 1]); colormap(ax(1),cbrewer('div','RdYlBu',100))
    axis off vis3d
    grid off
    view([1 0 0]);
    axis manual
    if slice_n==1; title('ODI');end
    
    ax(2)=subplot(length(apslices),4,2+(4*(slice_n-1))); cla
    hold on
    h = slice(xq,...
        yq,...
        zq,...
        vq_APsampled_space_azimuth,... % vq_inLGN_azimuth,... % 
        [apslices(slice_n)],...
        [],...
        []); % is plotted in terms of ARS for ease of navigation so RSA to ARS: z,x,
    alpha('color')
    set(h,'EdgeColor','none','FaceAlpha',1);
    set(gcf,'color','w')
        % plot dLGN outline
        plot3(ones(size(dLGN_boundary{slice_n})).*apslices(slice_n),slice_shape_points{slice_n}(dLGN_boundary{slice_n},2),...
            slice_shape_points{slice_n}(dLGN_boundary{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_ipsi{slice_n})).*apslices(slice_n),slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},2),...
            slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_shell{slice_n})).*apslices(slice_n),slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},2),...
            slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},3),'k','LineWidth',2)
        axis equal
        xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');
    ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right');
    zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
    caxis([0 90]); 
    colormap(ax(2),'jet');
%     colormap(ax(2),buildcmap('ygb'));
    
    axis off vis3d
    grid off
    view([1 0 0]);
    axis manual
    if slice_n==1; title('Azimuth');end
    
    ax(3)=subplot(length(apslices),4,3+(4*(slice_n-1))); cla
    hold on
    h = slice(xq,...
        yq,...
        zq,...
        vq_APsampled_space_elevation,... % vq_inLGN_elevation,... % 
        [apslices(slice_n)],...
        [],...
        []); % is plotted in terms of ARS for ease of navigation so RSA to ARS: z,x,
    alpha('color')
    set(h,'EdgeColor','none','FaceAlpha',1);
    set(gcf,'color','w')
        % plot dLGN outline
        plot3(ones(size(dLGN_boundary{slice_n})).*apslices(slice_n),slice_shape_points{slice_n}(dLGN_boundary{slice_n},2),...
            slice_shape_points{slice_n}(dLGN_boundary{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_ipsi{slice_n})).*apslices(slice_n),slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},2),...
            slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_shell{slice_n})).*apslices(slice_n),slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},2),...
            slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},3),'k','LineWidth',2)
    axis equal
    xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');
    ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right');
    zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
    caxis([-20 40]);
    colormap(ax(3),'jet')
%     colormap(ax(3),buildcmap('mc'));
    axis off vis3d
    grid off
    view([1 0 0]);
    axis manual
    if slice_n==1; title('Elevation');end
    
    ax(4)=subplot(length(apslices),4,4+(4*(slice_n-1))); cla
    hold on
    h = slice(xq,...
        yq,...
        zq,...
        idx_retmap_inbinoarea,... % vq_inLGN_azimuth,... % 
        [apslices(slice_n)],...
        [],...
        []); % is plotted in terms of ARS for ease of navigation so RSA to ARS: z,x,
        alpha('color')
        set(h,'EdgeColor','none','FaceAlpha',1);
        set(gcf,'color','w')
        % plot dLGN outline
        plot3(ones(size(dLGN_boundary{slice_n})).*apslices(slice_n),slice_shape_points{slice_n}(dLGN_boundary{slice_n},2),...
            slice_shape_points{slice_n}(dLGN_boundary{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_ipsi{slice_n})).*apslices(slice_n),slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},2),...
            slice_shape_points_ipsi{slice_n}(dLGN_boundary_ipsi{slice_n},3),'k','LineWidth',2)
        % plot dLGN ipsi outline
        plot3(ones(size(dLGN_boundary_shell{slice_n})).*apslices(slice_n),slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},2),...
            slice_shape_points_shell{slice_n}(dLGN_boundary_shell{slice_n},3),'k','LineWidth',2)
        axis equal
        xlim([min(shape_dLGN.Points(:,1)),max(shape_dLGN.Points(:,1))]); set(gca,'xDir','normal'); xlabel('x: anterior-posterior');
    ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right');
    zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
    caxis([0 1]); 
    colormap(ax(4),flip(cbrewer('seq','YlGn',100)));
    
    axis off vis3d
    grid off
    view([1 0 0]);
    axis manual
    if slice_n==1; title('Binocular area');end
        
    
    Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget','XLim','YLim','ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    
    % change to slice view
    view([1 0 0]); pause(0.1)
    xlim([apslices(displayed_slice)-2 apslices(displayed_slice)+2]);
    ylim([750,max(shape_dLGN.Points(:,2))]); set(gca,'yDir','normal'); ylabel('y: left-right');
    zlim([min(shape_dLGN.Points(:,3)),max(shape_dLGN.Points(:,3))]); set(gca,'zDir','reverse'); zlabel('z: superior-inferior');
    
   
end

for i = 1:4
    oldpos = ax(i).Position;
    cbar = colorbar(ax(i));
    cbar.Location = 'southoutside';
    cbar.FontSize = 10;
    ax(i).Position=oldpos;
end

subplot(ax(3))
plot3([max(xlim) max(xlim)],[max(ylim)-30 max(ylim)],[min(zlim) min(zlim)],'k','LineWidth',3)

%% plot interp eye preference map in the 2D retino space
nans_temp=isnan(vq_APsampled_space_odi(:)) |...
    isnan(vq_APsampled_space_elevation(:));
odi_temp = vq_APsampled_space_odi(:);
elevation_temp = vq_APsampled_space_elevation(:);
azimuth_tmep = vq_APsampled_space_azimuth(:);

odi_temp(nans_temp)=[];
elevation_temp(nans_temp)=[];
azimuth_tmep(nans_temp)=[];

[N,XEDGES,YEDGES,BINX,BINY]=histcounts2(elevation_temp,azimuth_tmep,'BinWidth',10);
acum_odi = accumarray([BINX,BINY],odi_temp,[],@(x) {x});
acum_odi_mean = cellfun(@(x) mean(x),acum_odi);
acum_odi_n = cellfun(@(x) length(x),acum_odi);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'eyepref_over_ret', 'Position', [1379,36,538,375]);
fig_number = fig_number+1;

h = imagesc(acum_odi_mean);
colormap(cbrewer('div','RdYlBu',100)); caxis([-1 1]);
set(h,'AlphaData',~isnan(acum_odi_mean))
set(gca,'YDir','normal')
xticks([0.5:size(acum_odi_mean,2)+0.5])
xticklabels(cellfun(@(x) num2str(x),num2cell(YEDGES),'UniformOutput', false))
yticks([0.5:size(acum_odi_mean,1)+0.5])
yticklabels(cellfun(@(x) num2str(x),num2cell(XEDGES),'UniformOutput', false))
axis equal
xlim([0.5 length(YEDGES)-0.5])
ylim([0.5 length(XEDGES)-0.5])
set(gca,'TickDir','out'); box off
cbar = colorbar;
cbar.Title.String = 'ODI';
xlabel('azimuth (°)')
ylabel('elevation (°)')
hold on

figure; set(gcf,'color','w','Renderer','opengl')
h = imagesc(acum_odi_mean);
colormap(cbrewer('div','RdYlBu',100)); caxis([-1 1]);
set(h,'AlphaData',~isnan(acum_odi_mean))
set(gca,'YDir','normal')
xticks([0.5:size(acum_odi_mean,2)+0.5])
xticklabels(cellfun(@(x) num2str(x),num2cell(YEDGES),'UniformOutput', false))
yticks([0.5:size(acum_odi_mean,1)+0.5])
yticklabels(cellfun(@(x) num2str(x),num2cell(XEDGES),'UniformOutput', false))
axis equal
xlim([0.5 length(YEDGES)-0.5])
ylim([0.5 length(XEDGES)-0.5])
set(gca,'TickDir','out'); box off
colorbar
xlabel('azimuth (°)')
ylabel('elevation (°)')
hold on
for i = 1:length(ccfpos_temp2)
    ccf_pos_temp = ccfpos_temp2(i,[3,1,2]);
    idx_temp = find(xq==ccf_pos_temp(1) & yq ==ccf_pos_temp(2) & zq==ccf_pos_temp(3));
    azim_pos(i)=vq_APsampled_space_azimuth(idx_temp);
    elev_pos(i)=vq_APsampled_space_elevation(idx_temp);
    cells_in_sampled_space_eyedom(i) = (ccfpos_temp2(i,4)>0)*2-1;
end
h = scatter(azim_pos/mean(diff(YEDGES))+find(YEDGES==0)-0.5,elev_pos/mean(diff(XEDGES))+find(XEDGES==0)-0.5,[],cells_in_sampled_space_eyedom,'filled','MarkerEdgeColor','k');
plot(binocular_boarder(:,1)/mean(diff(YEDGES))+find(YEDGES==0)-0.5,binocular_boarder(:,2)/mean(diff(XEDGES))+find(XEDGES==0)-0.5,'--k');

%% plot interp binocularity map in the 2D retino space
nans_temp=isnan(vq_APsampled_space_absodi(:)) |...
    isnan(vq_APsampled_space_elevation(:));
absodi_temp = vq_APsampled_space_absodi(:);
elevation_temp = vq_APsampled_space_elevation(:);
azimuth_tmep = vq_APsampled_space_azimuth(:);

absodi_temp(nans_temp)=[];
elevation_temp(nans_temp)=[];
azimuth_tmep(nans_temp)=[];

[N,XEDGES,YEDGES,BINX,BINY]=histcounts2(elevation_temp,azimuth_tmep,'BinWidth',10);
acum_absodi = accumarray([BINX,BINY],absodi_temp,[],@(x) {x});
acum_absodi_mean = cellfun(@(x) mean(x),acum_absodi);
acum_absodi_n = cellfun(@(x) length(x),acum_absodi);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'absODI_over_ret', 'Position', [1379,36,538,375]);
fig_number = fig_number+1;

h = imagesc(acum_absodi_mean);
colormap(flip(cbrewer('seq','YlGn',100))); caxis([0 1]);
set(h,'AlphaData',~isnan(acum_absodi_mean))
set(gca,'YDir','normal')
xticks([0.5:size(acum_absodi_mean,2)+0.5])
xticklabels(cellfun(@(x) num2str(x),num2cell(YEDGES),'UniformOutput', false))
yticks([0.5:size(acum_absodi_mean,1)+0.5])
yticklabels(cellfun(@(x) num2str(x),num2cell(XEDGES),'UniformOutput', false))
axis equal
xlim([0.5 length(YEDGES)-0.5])
ylim([0.5 length(XEDGES)-0.5])
set(gca,'TickDir','out'); box off
cbar = colorbar;
cbar.Title.String = '1-|ODI|';
xlabel('azimuth (°)')
ylabel('elevation (°)')
hold on


figure; set(gcf,'color','w','Renderer','opengl')
h = imagesc(acum_absodi_mean);
colormap(flip(cbrewer('seq','YlGn',100))); caxis([0 1]);
set(h,'AlphaData',~isnan(acum_absodi_mean))
set(gca,'YDir','normal')
xticks([0.5:size(acum_absodi_mean,2)+0.5])
xticklabels(cellfun(@(x) num2str(x),num2cell(YEDGES),'UniformOutput', false))
yticks([0.5:size(acum_absodi_mean,1)+0.5])
yticklabels(cellfun(@(x) num2str(x),num2cell(XEDGES),'UniformOutput', false))
axis equal
xlim([0.5 length(YEDGES)-0.5])
ylim([0.5 length(XEDGES)-0.5])
set(gca,'TickDir','out'); box off
cbar = colorbar;
cbar.Title.String = '1-|ODI|';
xlabel('azimuth (°)')
ylabel('elevation (°)')
hold on
for i = 1:length(ccfpos_temp2)
    ccf_pos_temp = ccfpos_temp2(i,[3,1,2]);
    idx_temp = find(xq==ccf_pos_temp(1) & yq ==ccf_pos_temp(2) & zq==ccf_pos_temp(3));
    azim_pos(i)=vq_APsampled_space_azimuth(idx_temp);
    elev_pos(i)=vq_APsampled_space_elevation(idx_temp);
    cells_in_sampled_space_ODI(i) = ccfpos_temp2(i,4);
    cells_in_sampled_space_absODI(i) = abs(ccfpos_temp2(i,4));
end
h = scatter(azim_pos/mean(diff(YEDGES))+find(YEDGES==0)-0.5,elev_pos/mean(diff(XEDGES))+find(XEDGES==0)-0.5,[],1-cells_in_sampled_space_absODI,'filled','MarkerEdgeColor','k');
plot(binocular_boarder(:,1)/mean(diff(YEDGES))+find(YEDGES==0)-0.5,binocular_boarder(:,2)/mean(diff(XEDGES))+find(XEDGES==0)-0.5,'--k');

%% ODI of neurons in the dLGN subserving the Dräger defined binocular visual field
idx_inbinoarea=inpolygon(azim_pos,elev_pos,...
    [binocular_boarder(:,1);flip(binocular_boarder(:,1)*-1)],...
    [binocular_boarder(:,2);flip(binocular_boarder(:,2))]);

temp_absODI = cells_in_sampled_space_absODI;
temp_ODI = cells_in_sampled_space_ODI;

figure; subplot(1,2,1); set(gcf, 'color','w')
scatter(azim_pos(idx_inbinoarea),elev_pos(idx_inbinoarea),[],'filled','MarkerFaceColor','m','MarkerEdgeColor','m'); hold on
scatter(azim_pos(~idx_inbinoarea),elev_pos(~idx_inbinoarea),[],'filled','MarkerFaceColor','g','MarkerEdgeColor','g');
plot([binocular_boarder(:,1);flip(binocular_boarder(:,1)*-1)],...
    [binocular_boarder(:,2);flip(binocular_boarder(:,2))],'--k')
axis equal tight
xlabel('azimuth (°)')
ylabel('elevation (°)')

subplot(1,2,2)
plotSpread({temp_absODI(idx_inbinoarea),temp_absODI(~idx_inbinoarea)},...
    'categoryMarkers',{'.','.'},...
    'distributionColors',{'m','g'}); hold on 

plot([1],median(temp_absODI(idx_inbinoarea)),'sr');
plot([2],median(temp_absODI(~idx_inbinoarea)),'sr');
plot([1],mean(temp_absODI(idx_inbinoarea)),'sk');
plot([2],mean(temp_absODI(~idx_inbinoarea)),'sk');

[p1,~,STATS1]=ranksum(temp_absODI(idx_inbinoarea),temp_absODI(~idx_inbinoarea),'tail', 'both');

title({['binoarea, mean: ' num2str(mean(temp_absODI(idx_inbinoarea)),2) ...
    ', median: ' num2str(median(temp_absODI(idx_inbinoarea)),2)];...
    ['monoarea, mean: ' num2str(mean(temp_absODI(~idx_inbinoarea)),2) ...
    ', median: ' num2str(median(temp_absODI(~idx_inbinoarea)),2)];...
    ['MWU (n=' num2str(sum(idx_inbinoarea)) ' & ' num2str(sum(~idx_inbinoarea))...
    ', p=' num2str(p1,2), ', U=' num2str(STATS1.ranksum) ')']});

xticklabels({['bino\newlinezone'],['mono\newlinezone'],['shell'],['core']})
set(gca,'TickLabelInterpreter','tex')
ylabel('|ODI|')

idx_inbinoarea_ns = find(idx_inbinoarea);
idx_inmonoarea_ns = find(~idx_inbinoarea);


fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'ODI_dist_in bino_area', 'Position', [1,37,506,287]); 
fig_number = fig_number+1;

clear ax
ax(1) = subplot(1,2,1); 
binos = abs(temp_ODI(idx_inbinoarea_ns))<1; 
binobincounts = histcounts(temp_ODI(idx_inbinoarea_ns(binos)),[-1:2/6:1]); 
binscounts = [sum(temp_ODI(idx_inbinoarea_ns)==-1), 0, binobincounts, 0, sum(temp_ODI(idx_inbinoarea_ns)==1)];
b = bar(binscounts/sum(binscounts)*100,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title({'binocular region';['n=' num2str(sum(binscounts))]})

ax(2) = subplot(1,2,2); 
binos = abs(temp_ODI(idx_inmonoarea_ns))<1; 
binobincounts = histcounts(temp_ODI(idx_inmonoarea_ns(binos)),[-1:2/6:1]); 
binscounts = [sum(temp_ODI(idx_inmonoarea_ns)==-1), 0, binobincounts, 0, sum(temp_ODI(idx_inmonoarea_ns)==1)];
b = bar(binscounts/sum(binscounts)*100,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title({'monocular region';['n=' num2str(sum(binscounts))]})

for i = 1:length(ax)
    set(ax(i),'XTick', [1,2.5,5.5,8.5 10])
    set(ax(i),'XTickLabels', {'-1','>-1','0','<1','1'})
    set(ax(i),'box','off')
    ax(i).YLabel.String='Fraction of cells %';
    set(ax(i),'TickDir','out')
    ax(i).XLabel.String='ODI';
end

%% shell vs core vs vlgncomparison
filter_dLGN = find(QC_cell_filter(data, your_data_dir,...
    'transduction_QC', 1, ...
    'dLGN', 1, ...
    'A_ramp_signal_QC',1,...
    'CS_internal',1,...
    'Interneuron_morph',0,...
    'CCF_QC',1));

ODI_A_dlgn_odi = [data(filter_dLGN).ODI_AMPA_step_peak];

dLGN_shell_odi = find([data(filter_dLGN).dLGN_shell]);
dLGN_core_odi = find(~[data(filter_dLGN).dLGN_shell]);
dLGN_contra_zone_odi = find([data(filter_dLGN).dLGN_contra_zone]);
dLGN_ipsi_zone_odi = find([data(filter_dLGN).dLGN_ipsi_zone]);

fig_hand(fig_number) = figure;
set(fig_hand(fig_number), 'Name', 'core_vs_shell_ODI', 'Position', [511,37,977,287]); 
fig_number = fig_number+1;

clear ax
ax(1) = subplot(1,4,1); 
binos = abs(ODI_A_dlgn_odi(dLGN_contra_zone_odi))<1; 
binobincounts = histcounts(ODI_A_dlgn_odi(dLGN_contra_zone_odi(binos)),[-1:2/6:1]); 
binscounts = [sum(ODI_A_dlgn_odi(dLGN_contra_zone_odi)==-1), 0, binobincounts, 0, sum(ODI_A_dlgn_odi(dLGN_contra_zone_odi)==1)];
b = bar(binscounts/sum(binscounts)*100,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title({'Contra core';['n=' num2str(sum(binscounts))]})

ax(2) = subplot(1,4,2); 
binos = abs(ODI_A_dlgn_odi(dLGN_ipsi_zone_odi))<1; 
binobincounts = histcounts(ODI_A_dlgn_odi(dLGN_ipsi_zone_odi(binos)),[-1:2/6:1]); 
binscounts = [sum(ODI_A_dlgn_odi(dLGN_ipsi_zone_odi)==-1), 0, binobincounts, 0, sum(ODI_A_dlgn_odi(dLGN_ipsi_zone_odi)==1)];
b = bar(binscounts/sum(binscounts)*100,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title({'Ipsi core';['n=' num2str(sum(binscounts))]})

ax(3) = subplot(1,4,3); 
binos = abs(ODI_A_dlgn_odi(dLGN_shell_odi))<1; 
binobincounts = histcounts(ODI_A_dlgn_odi(dLGN_shell_odi(binos)),[-1:2/6:1]); 
binscounts = [sum(ODI_A_dlgn_odi(dLGN_shell_odi)==-1), 0, binobincounts, 0, sum(ODI_A_dlgn_odi(dLGN_shell_odi)==1)];
b = bar(binscounts/sum(binscounts)*100,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title({'Shell';['n=' num2str(sum(binscounts))]})

ax(4) = subplot(1,4,4); 
binos = abs(ODI_A_dlgn_odi(dLGN_core_odi))<1; 
binobincounts = histcounts(ODI_A_dlgn_odi(dLGN_core_odi(binos)),[-1:2/6:1]); 
binscounts = [sum(ODI_A_dlgn_odi(dLGN_core_odi)==-1), 0, binobincounts, 0, sum(ODI_A_dlgn_odi(dLGN_core_odi)==1)];
b = bar(binscounts/sum(binscounts)*100,1);
b.FaceColor = 'flat';
b.CData(:,:)=[colorscale_ODI(1,:); [0 0 0]; colorscale_ODI(2:2:end,:); [0 0 0]; colorscale_ODI(end,:)];
title({'Core';['n=' num2str(sum(binscounts))]})

for i = 1:length(ax)
    set(ax(i),'XTick', [1,2.5,5.5,8.5 10])
    set(ax(i),'XTickLabels', {'-1','>-1','0','<1','1'})
    set(ax(i),'box','off')
    ax(i).YLabel.String='Fraction of cells %';
    set(ax(i),'TickDir','out')
    ax(i).XLabel.String='ODI';
end

[tbl,chi2,p] = crosstab([ODI_A_dlgn_odi(dLGN_contra_zone_odi)>0,ODI_A_dlgn_odi(dLGN_ipsi_zone_odi)>0],...
    [ones(1,length(dLGN_contra_zone_odi)),zeros(1,length(dLGN_ipsi_zone_odi))]+2);

array2table(tbl, 'RowNames',{'ipsi_dom' 'contra_dom'}, 'VariableNames',{'ipsi_zone' 'contra_zone'} )
p
disp(['missplaced cells: ' num2str((tbl(1,2)+tbl(2,1))/sum(tbl(:)))])

figure;
plotSpread({1-abs(ODI_A_dlgn_odi(dLGN_contra_zone_odi)),...
    1-abs(ODI_A_dlgn_odi(dLGN_ipsi_zone_odi)),...
    1-abs(ODI_A_dlgn_odi(dLGN_shell_odi)),...
    1-abs(ODI_A_dlgn_odi(dLGN_core_odi))},...
    'categoryMarkers',{'.','.','.','.'},...
    'distributionColors',{'k','k','k','k'}); hold on 
plot([1],median(1-abs(ODI_A_dlgn_odi(dLGN_contra_zone_odi))),'sr');
plot([2],median(1-abs(ODI_A_dlgn_odi(dLGN_ipsi_zone_odi))),'sr');
plot([3],median(1-abs(ODI_A_dlgn_odi(dLGN_shell_odi))),'sr');
plot([4],median(1-abs(ODI_A_dlgn_odi(dLGN_core_odi))),'sr');
xticklabels({['contra\newlinezone'],['ipsi\newlinezone'],['shell'],['core']})
set(gca,'TickLabelInterpreter','tex')
ylabel('1-|ODI|')

[p1,~,STATS1]=ranksum(1-abs(ODI_A_dlgn_odi(dLGN_contra_zone_odi)),1-abs(ODI_A_dlgn_odi(dLGN_ipsi_zone_odi)),'tail', 'both');
[p2,~,STATS2]=ranksum(1-abs(ODI_A_dlgn_odi(dLGN_shell_odi)),1-abs(ODI_A_dlgn_odi(dLGN_core_odi)),'tail', 'both');

title({['MWU (n=' num2str(length(dLGN_contra_zone_odi)) ' & ' num2str(length(dLGN_ipsi_zone_odi))...
    ', p=' num2str(p1,2), ', U=' num2str(STATS1.ranksum) ')'],...
    ['MWU (n=' num2str(length(dLGN_shell_odi)) ' & ' num2str(length(dLGN_core_odi))...
    ', p=' num2str(p2,2) ', U=' num2str(STATS2.ranksum) ')']})

%% save figs
if savefigs
    mkdir([savefigs_location '\allenfigs\'])
    for i=1:length(fig_hand)
        saveas(fig_hand(i),[savefigs_location '\allenfigs\' fig_hand(i).Name],'svg')
        if i == 1 
            temp = findall(fig_hand(i),'type','axes');
            for ii = 1:length(temp)
                temp(ii).Visible = 'off';
            end
            saveas(fig_hand(i),[savefigs_location '\allenfigs\' fig_hand(i).Name],'png')
            for ii = 1:length(temp)
                temp(ii).Visible = 'on';
            end
        end
        saveas(fig_hand(i),[savefigs_location '\allenfigs\' fig_hand(i).Name],'fig');
    end
end