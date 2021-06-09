% produce and save final figs
%% load main data structure
your_data_dir = 'I:\Martin Fernholz\LGN _Project_Common\code and data upload\data\';
cd(your_data_dir)
d = dir('Full_Data_*.mat');
file_name = d(find([d(:).datenum]==max([d(:).datenum]))).name;
data = load_StrArray(file_name,cd,'IncludeField','all');

%% load ACCF data and retinotopy
load([your_data_dir 'Allen_and_retino_data.mat'])

%% set fig defaults
set(groot,'DefaultAxesTickDir', 'out')
set(groot,'DefaultAxesTickDirMode', 'manual');
set(groot,'DefaultAxesFontWeight', 'normal');
set(groot,'DefaultAxesFontName', 'Arial');
set(groot,'DefaultAxesFontSizeMode', 'manual')
set(groot,'DefaultAxesFontSize', 10)
set(groot,'DefaultAxesTitleFontSizeMultiplier',1.2)
set(groot,'DefaultAxesTickLength', [0.02 0.025]);
set(groot,'DefaultFigureColor','w')
set(groot,'DefaultLineLineWidth',2)
set(groot,'DefaultHistogramFaceColor',[0.5 0.5 0.5])
set(groot,'DefaultScatterMarkerEdgeColor',[0.5 0.5 0.5])
set(groot,'DefaultScatterMarkerFaceColor',[0.5 0.5 0.5])
set(groot,'DefaultPolaraxesFontSizeMode', 'manual')
set(groot,'DefaultPolaraxesFontWeight', 'normal');
set(groot,'DefaultAxesFontName', 'Arial');
set(groot,'DefaultPolaraxesFontSize',10)
set(groot,'DefaultPolaraxesTitleFontSizeMultiplier',1.2);
set(groot,'DefaultTextFontName', 'Arial');
set(groot,'DefaultLegendFontName', 'Arial');
set(groot,'DefaultColorbarFontName', 'Arial');
set(groot,'DefaultColorbarTickDirection','out');
set(groot,'DefaultFigureRenderer','painter');

savefigs=0;
savefigs_location = []; % Enter destination

%% create plots
% fig1&2 Overview
categorical_vs_continuous_binoc

% fig3 + supfigs: allen ccf and retinotopy
allenplots_figs

% fig6 + supfigs: morph overview and effect on monocularity
morph_figs
morph_clustering_figs
dend_bias_figs

% fig7 + subfigs: synaptic input selection and refinement
FR_intro_figs
synaptic_mech_figs
