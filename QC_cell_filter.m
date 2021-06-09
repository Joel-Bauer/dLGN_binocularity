function filter_out = QC_cell_filter(data, QC_loc,  varargin)
% this will be the function were you choose what criteria cells need to
% fullfill to be included. inputs will be name value pairs of condisitons
% (eg. 'brain_contra_ipsi',1,'moropho_QC',1).
% varargin = {'brain_contra_ipsi',0};
% possible filters and their default values:
% filters based on data structure (and default values)
%     brain_contra_ipsi: nan % when the filter is set to nan it is not used 
%     dLGN: nan
%     dLGN_shell: nan
%     dLGN_contra_zone: nan
%     dLGN_ipsi_zone: nan
%     vLGN: nan
% filters based on QC sheet (and default values)
%     transduction_QC: 1 % indicate the min rating of transduction 
%             0=not good enough
%             1=good local transduction 
%             2=ok overall transduction
%             3=perfect transduction
%     Conf_QC: nan
%     morph_QC: nan
%     TwoPTrans_QC: nan
%     CCF_QC: nan
%     A_ramp_signal_QC: nan
%     N_ramp_signal_QC: nan
%     AtoN_Rschange_QC: nan
%     CS_internal: nan % if 0 then also checks for CC data
%     Interneuron_morph: nan

if isempty(QC_loc)
    QC_loc = 'I:\Martin Fernholz\LGN _Project_Common\code and data upload\data\';
end

p = inputParser();
p.KeepUnmatched = false;
p.CaseSensitive = false;
p.StructExpand  = false;
p.PartialMatching = true;

% filters based on data structure
addParameter(p, 'brain_contra_ipsi' , nan) % when the filter is set to nan it is not used 
addParameter(p, 'dLGN', nan); 
addParameter(p, 'dLGN_shell', nan);
addParameter(p, 'dLGN_contra_zone', nan);
addParameter(p, 'dLGN_ipsi_zone', nan);
addParameter(p, 'vLGN', nan);

% filters based on QC sheet
addParameter(p, 'transduction_QC', 1); % 
addParameter(p, 'Conf_QC', nan);
addParameter(p, 'morph_QC', nan);
addParameter(p, 'TwoPTrans_QC', nan);
addParameter(p, 'CCF_QC', nan);
addParameter(p, 'A_ramp_signal_QC', nan);
addParameter(p, 'N_ramp_signal_QC', nan);
addParameter(p, 'AtoN_Rschange_QC', nan);
addParameter(p, 'CS_internal', nan);
addParameter(p, 'Interneuron_morph', 0); % should be 0

parse(p, varargin{:})
cell_filters = p.Results;

% remove fields with nan
fn = fieldnames(cell_filters);
for k=1:numel(fn)
    if isnumeric(cell_filters.(fn{k})) && isnan(cell_filters.(fn{k}))
        cell_filters = rmfield(cell_filters,fn{k});
    end
end

% get QC sheet
originaldir = cd;
cd(QC_loc)
d = dir('QC_Filter*.xlsx');
file_name = d(find([d(:).datenum]==max([d(:).datenum]))).name;
QC_table = readtable([file_name]);
QC_table = table2struct(QC_table);
cd(originaldir) % go back to where you came from!

data_cell_names = join([cellfun(@(x) x(1:6),{data(:).patching_date},'UniformOutput',false)' ...
    {data(:).Setup}'...
    {data(:).cellname}']);
QC_table_cell_names = join([{QC_table(:).file_name}' ...
    {QC_table(:).Setup}'...
    {QC_table(:).cellname}']);

% remove cells not present in data from QC_table
try
    QC_table(~ismember(unique(flip(QC_table_cell_names),unique(data_cell_names)))) = [];
end

% match order of QC_table and add QC subfields to data (do not overwrite those already present!)
I = cell2mat(cellfun(@(x) strmatch(x,data_cell_names),QC_table_cell_names, 'UniformOutput', false));
fn_QC = fieldnames(QC_table);
fn_data = fieldnames(data);
data2 = data; 
for fn_idx = 6:length(fn_QC)
    for celln_QC = 1:length(I)
        celln_data = I(celln_QC);
        data2(celln_data).(fn_QC{fn_idx}) = QC_table(celln_QC).(fn_QC{fn_idx});
    end
    for celln_data = 1:length(data)
        if isempty(data2(celln_data).(fn_QC{fn_idx}))
            data2(celln_data).(fn_QC{fn_idx}) = nan;
        end
    end
end

% loop through filter and create filter vectors
fn = fieldnames(cell_filters);
filters_array = [];
for k = 1:length(fn)
    if strcmp(fn{k},'dLGN')
        temp_filter = cell_filters.dLGN == [([data2.dLGN_contra_zone] + [data2.dLGN_ipsi_zone])==1];
        filters_array(k,:) = temp_filter;
        
    elseif strcmp(fn{k},'transduction_QC')
        temp_filter = cell_filters.transduction_QC <= [data2.transduction_QC];
        filters_array(k,:) = temp_filter;
        
    elseif strcmp(fn{k},'Conf_QC')
        temp_filter = cellfun(@(x) all(cell_filters.Conf_QC==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
        
    elseif strcmp(fn{k},'TwoPTrans_QC')
        temp_filter = cellfun(@(x) all(cell_filters.TwoPTrans_QC==x),{data2.(fn{k})});
        temp_filter2 = cellfun(@(x) ~isempty(x.rotation), {data2.TwoP_to_conf_transform});
        
        if cell_filters.TwoPTrans_QC == 1
            filters_array(k,:) = temp_filter & temp_filter2;
        elseif cell_filters.TwoPTrans_QC == 0
            filters_array(k,:) = temp_filter | ~temp_filter2;
        end
        
    elseif strcmp(fn{k},'CCF_QC')
        temp_filter = cellfun(@(x) all(cell_filters.CCF_QC==x),{data2.(fn{k})});
        temp_filter2 = cellfun(@(x) ~isempty(x), {data2.ccfv3_pos});
        
        if cell_filters.CCF_QC == 1
            filters_array(k,:) = temp_filter & temp_filter2;
        elseif cell_filters.CCF_QC == 0
            filters_array(k,:) = temp_filter | ~temp_filter2;
        end
        
    elseif strcmp(fn{k},'CS_internal')
        temp_filter = cellfun(@(x) all(cell_filters.CS_internal==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
        
    elseif strcmp(fn{k},'morph_QC')
        min_total_dend_length = 600;
        temp_filter = cellfun(@(x) all(cell_filters.morph_QC==x),{data2.(fn{k})});
%         temp_filter2 = cellfun(@(x) ~isempty(x), {data2.interp_morphology});
        temp_filter2 = cellfun(@(x) ~isempty(x.interp_morphology), {data2.morphology});
        temp_filter3= ones(1,length(temp_filter2));
        temp = find(temp_filter2);
        temp2 = cellfun(@(x) x.Sholl_analysis.total_length<min_total_dend_length, {data2(temp).morphology});
        temp = temp(temp2);
        temp_filter3(temp) = 0;
        
        if cell_filters.morph_QC == 1
            filters_array(k,:) = temp_filter & temp_filter2 & temp_filter3;
        elseif cell_filters.morph_QC == 0
            filters_array(k,:) = temp_filter | ~temp_filter2 | ~temp_filter3;
        end
        
    elseif strcmp(fn{k},'A_ramp_signal_QC')
        temp_filter = cellfun(@(x) all(cell_filters.A_ramp_signal_QC==x),{data2.(fn{k})});
        temp_filter2 = cellfun(@(x) ~isempty(x), {data2.ODI_AMPA_step_peak});
        
        if cell_filters.A_ramp_signal_QC == 1
            filters_array(k,:) = temp_filter & temp_filter2;
        elseif cell_filters.A_ramp_signal_QC == 0
            filters_array(k,:) = temp_filter | ~temp_filter2;
        end
        
    elseif strcmp(fn{k},'N_ramp_signal_QC')
        temp_filter = cellfun(@(x) all(cell_filters.N_ramp_signal_QC==x),{data2.(fn{k})});
        temp_filter2 = cellfun(@(x) ~isempty(x), {data2.ODI_NMDA_step_peak});
        
        if cell_filters.N_ramp_signal_QC == 1
            filters_array(k,:) = temp_filter & temp_filter2;
        elseif cell_filters.N_ramp_signal_QC == 0
            filters_array(k,:) = temp_filter | ~temp_filter2;
        end
        
    else
        temp_filter = cellfun(@(x) all(cell_filters.(fn{k})==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
        
    end
end
if isempty(filters_array)
    filter_out = ones(1,length(data2));
else
    filter_out = all(filters_array,1);
end
disp([num2str(sum(filter_out)) ' cells pass all filters' ])
    
    
    