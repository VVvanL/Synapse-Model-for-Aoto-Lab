% pilot script to convert model output (open channel number) to mock IPSCs
clearvars; close all

%% ==========================================
Vm = -70; % resting membrane potential (in mV)
Vr = -85; % reversal potential (in mV)
g = 0.027; % single-channel conductance (in msec)
%% ==========================================

folderN = uigetdir; folderN = [folderN,filesep];
foldparts = strsplit(folderN,filesep); dirname = foldparts{end-1}; clear foldparts

IPSC = struct(); % create data structure for storing new values

% get list of files, assumes only model output mat files are in directory
filelist = dir(folderN); filelist(1:2) = []; file_n = size(filelist,1); 

for f = 1:file_n
    rootname = filelist(f).name(1:end-4); % file name without extension
    load([folderN, filelist(f).name]);

    n = data_matrix.R_open; % extract matrix of open receptors
    % n = open_smoothed;
    I = (g * n) * (Vm - Vr); % matrix operations to retrieve mock IPSC curve
    % save new variable to file
    save([folderN, filelist(f).name], 'I', '-append')
    % add IPSC data to structure
    IPSC.(rootname).I = I;

    %TODO: Add fitting of rise and decay    
end



params = struct();
params.Vm = Vm;
params.Vr = Vr;
params.g = g;
params.time_sec = time_sec;

save([folderN, dirname, '_IPSC_data.mat'], 'params', 'IPSC' )

model_names = fieldnames(IPSC); model_n = length(model_names);

%% create inverse of IPCS (convert to negative values)
for md = 1:model_n
    mdl_str = model_names{md};    
    IPSC.(mdl_str).I_neg = -IPSC.(mdl_str).I;
end
save([folderN, dirname, '_IPSC_data.mat'], 'params', 'IPSC' )

%% create fractional (proportion) normalized IPSC curve for each type
% type fraction in each compartment (create vector to match model names order)
type_fraction = [0.12, 0.88, 0.42, 0.58];
for md = 1:model_n
    mdl_str = model_names{md};
    params.(mdl_str).fraction = type_fraction(md);
    IPSC.(mdl_str).fractional_norm = IPSC.(mdl_str).I .* type_fraction(md);
end

% proportionately convolve compartment models
compartment_str = {'dendritic', 'somatic'}; cmp_n = length(compartment_str);
model_cmp = [1,2;3,4];  % models corresponding to compartment

for cmp = 1:cmp_n
    cmp_str = compartment_str{cmp};
    model_numbers = model_cmp(cmp,:);
    IPSC.(cmp_str).fractionally_convolved = ...
        IPSC.(model_names{model_numbers(1)}).fractional_norm + IPSC.(model_names{model_numbers(2)}).fractional_norm;

end


%% Plotting functions (optional call)
% TODO: replace shadedErrorBar function with more flexible plotting scheme
figure;
t = tiledlayout(1,2); 
for f = 1:model_n
    title_str = model_names{f};
    nexttile; hold on; grid on

    shadedErrorBar(params.time_sec,mean(IPSC.(title_str).I_neg, 2),std(IPSC.(title_str).I_neg, 0, 2),...
        'lineprops','-k','transparent',true);
    ylim([-17 0]); 
    title(title_str, 'Interpreter','none')

end

% plot fractionally convolved IPSCs per compartment\
figure;
t2 = tiledlayout(1,2);
for cmp = 1:cmp_n
    title_str = compartment_str{cmp};
    nexttile; hold on; grid on

    shadedErrorBar(params.time_sec,mean(IPSC.(title_str).fractionally_convolved, 2),...
        std(IPSC.(title_str).fractionally_convolved, 0, 2),'lineprops','-k','transparent',true);
    ylim([0 25]); 
    title(title_str, 'Interpreter','none')

end



% plot individual fraction norms
figure;
t3 = tiledlayout(2,2); 
for f = 1:model_n
    title_str = model_names{f};
    nexttile; hold on; grid on

    shadedErrorBar(params.time_sec,mean(IPSC.(title_str).fractional_norm, 2),std(IPSC.(title_str).fractional_norm, 0, 2),...
        'lineprops','-k','transparent',true);
    ylim([0 14]); 
    title(title_str, 'Interpreter','none')

end
