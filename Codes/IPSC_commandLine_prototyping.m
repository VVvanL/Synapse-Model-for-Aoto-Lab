% File containing command line (not script ready) prototyping for IPSC fitting

%% get file name/directory information
[file, path] = uigetfile();
load([path,file])

model_names = fieldnames(IPSC); model_n = length(model_names);

% create mean trace for fitting
for md = 1:model_n
    mdl_str = model_names{md};    
    
    IPSC.(mdl_str).I_mean = mean(IPSC.(mdl_str).I, 2);

    ts_n = length(IPSC.(mdl_str).I_mean);
    I_max = max(IPSC.(mdl_str).I_mean);
    max_idx = find(IPSC.(mdl_str).I_mean == I_max);

    rise_phase = 1:max_idx; rise_phase = rise_phase';
    decay_phase = max_idx:ts_n; decay_phase = decay_phase';

    rise_phase = IPSC.(mdl_str).I_mean(rise_phase,:);

end
save([path, file], 'IPSC', '-append')

figure;
grid on; hold on
plot( IPSC.(mdl_str).I_mean,'k-','LineWidth', 1)
xlim([0 10e4])
plot(max_idx, I_max, 'm*')

figure; hold on; grid on
plot(rise_phase,'k','LineWidth', 1)

%% plot comparative rise phases
figure; hold on
for md = 1:model_n
    mdl_str = model_names{md};
    
    I_max = max(IPSC.(mdl_str).I_mean);
    max_idx = find(IPSC.(mdl_str).I_mean == I_max);

    rise_phase = 1:max_idx; rise_phase = rise_phase';
    decay_phase = max_idx:ts_n; decay_phase = decay_phase';

    rise_phase = IPSC.(mdl_str).I_mean(rise_phase,:);

    plot(rise_phase, 'LineWidth', 1.5)


end
legend(model_names, 'Interpreter','none')

figure; hold on
for md = 1:model_n
    mdl_str = model_names{md};
    
    I_max = max(IPSC.(mdl_str).I_mean);
    max_idx = find(IPSC.(mdl_str).I_mean == I_max);

    rise_phase = 1:max_idx; rise_phase = rise_phase';
    decay_phase = max_idx:ts_n; decay_phase = decay_phase';

    decay_phase = IPSC.(mdl_str).I_mean(decay_phase,:);

    plot(decay_phase, 'LineWidth', 1.5)


end
legend(model_names, 'Interpreter','none')