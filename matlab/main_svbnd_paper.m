
% include dependencies
addpath(genpath('integrators'));
addpath(genpath('problems'));

epsilons  = [1];    % diffusion coefficient
ic_types  = [2, 3]; % initial conditions (2 gaussian, 3 gaussian)
orders    = [2, 3]; % method orders
Nts       = round(logspace(1.3, 5, 20)); % num timesteps
use_cache = true; % set to false to re-run experiments

for i = 1 : length(epsilons)
    for j = 1 : length(ic_types)        
        for k = 1 : length(orders)        
            runExperiment(epsilons(i), ic_types(j), orders(k), Nts, use_cache);
        end
    end
end

function runExperiment(epsilon, ic_type, order, Nts, use_cache)

    % --> init prb and integrators
    prb = initProblem(epsilon, ic_type);
    [integrators, line_styles, ms, leg_ord] = initIntegrators(order);
    
    % --> set save path    
    base_dir  = 'figures/';
    base_name = replace(sprintf('nld-burgers-eps-%.2f-ic-%i-ord-%i', epsilon, ic_type, order),'.','_');
    
    % --> run experiment
    results_path = [base_dir, 'results-', base_name, '.mat'];
    if(use_cache && isfile(results_path)) % -> load existing results
        load(results_path, 'error', 'time')
    else % --> run experiment and save results
        [error, ~, time] = convergence(integrators, prb, Nts);
        save([base_dir, 'results-', base_name, '.mat'], 'prb', 'integrators', 'line_styles', 'error', 'time');
    end
    
    % re-order methods
    [integrators, line_styles, error, time] = reorderResults(integrators, line_styles, error, time, leg_ord);

    % --> generate figures
    base_path = [base_dir, base_name];

    % ----> convergence
    f = makeConvergencePlot(prb.tspan, Nts, error, line_styles, []);
    exportFigure(f, struct('SavePath', [base_path, '-convergence'], 'PaperPosition', 2 * [0 0 8 5], 'Format', 'pdf', 'Renderer', 'painters'));

    % ----> efficiency
    reg_str = '^ssp\d-([2-9]|(\d\d+)x)'; % remove composed explicit ssp schemes
    sub_inds_eff = cellfun(@(i) isempty(regexp(i.name,reg_str, 'once')), integrators);
    f = makeEfficiencyPlot(prb.tspan, time(:,sub_inds_eff), error(:,sub_inds_eff), [], line_styles(sub_inds_eff));
    exportFigure(f, struct('SavePath', [base_path, '-efficiency'], 'PaperPosition', 2 * [0 0 8 5], 'Format', 'pdf', 'Renderer', 'painters'));
    
    % ----> legend
    leg_str = cellfun(@(m) m.name, integrators, 'UniformOutput', false);
    f = makeLegendPlot(line_styles, leg_str);
    exportFigure(f, struct('SavePath', [base_path, '-legend'], 'PaperPosition', 2 * [0 0 2*8 5], 'Format', 'pdf', 'Renderer', 'painters'));
    
    % ----> solution
    base_name_sol = replace(sprintf('nld-burgers-eps-%.2f-ic-%i', epsilon, ic_type),'.','_');
    base_path_sol = [base_dir, base_name_sol];
    f = makeSolutionPlot(prb);
    exportFigure(f, struct('SavePath', [base_path_sol, '-solution'], 'PaperPosition', 2 * [0 0 8 5], 'Format', 'pdf', 'Renderer', 'opengl'));

end

function prb = initProblem(epsilon, ic_type)

    prb = SVBND();
    prb.epsilon = epsilon;
    prb.initial_condition_type = ic_type;
    prb.reset();

end

function [integrators, line_styles, ms, plot_order] = initIntegrators(order)

    switch(order)
        case 2
            ms = [1 4 16];
            [integrators, line_styles, plot_order] = methodsOrder2(ms);
        case 3
            ms = [1 2 4 8];
            [integrators, line_styles, plot_order] = methodsOrder3(ms, 1);
    end

end

function [methods, markers, legend_order] = methodsOrder3(ms, variant)
% Methods to run for third-order experiment

    if(nargin < 2)
        variant = 1;
    end

    mrnprk_constuctor = @(m) MR_NPRK_IMEX_SDIRK3_WRAP(@ssp33, m, variant, 2, sprintf('MR-NPRK3-%i[ssp3-%ix]', variant, m));
    erk_constructor   = @(m) ERK(@() nStepTableau(@ssp33, m), sprintf('ssp3-%ix', m));
    
    methods_nprk   = { NPRK(@NP_IMEX_354_Sa) };
    methods_mrnprk = arrayfun(mrnprk_constuctor, ms, 'UniformOutput', false);
    methods_erk    = arrayfun(erk_constructor, ms, 'UniformOutput', false);
    methods        = [methods_nprk, methods_mrnprk, methods_erk];

    cb             = .2 * [1 1 1]; % grey color for nprk (non multirate)
    dc             = defaultColors();
    
    markers_nprk   = {ls('square-.', dc(5,:), 5)};
    markers_mrnprk = {ls('<-', dc(1,:)), ls('>-', dc(2,:)), ls('*-', dc(3,:), 8), ls('diamond-', dc(4,:)), ls('o-', dc(5,:), 6), ls('^-', dc(6,:)), ls('hexagram-', dc(7,:))};
    markers_erk    = {ls('<:', cb),      ls('>:', cb),      ls('*:', cb, 8),      ls('diamond:', cb),      ls('o:', cb, 6),      ls('^:',cb),       ls('hexagram:', cb)};
    
    nm      = length(ms);
    markers = [markers_nprk, markers_mrnprk(1:nm), markers_erk(1:nm)];

    legend_order = [2 : (2 * nm + 1), 1];
    
end

function [methods, markers, legend_order] = methodsOrder2(ms)
% Methods to run for second-order experiment

    mrnprk_constuctor = @(m) MR_NPRK_IMEX_SDIRK2_WRAP(@ssp22, m, -1, sprintf('MR-NPRK2m-[ssp2-%ix]', m));
    erk_constructor   = @(m) ERK(@() nStepTableau(@ssp22, m), sprintf('ssp2-%ix', m));
    
    methods_nprk   = { NPRK(@NP_IMEX_242b) };
    methods_mrnprk = arrayfun(mrnprk_constuctor, ms, 'UniformOutput', false);
    methods_erk    = arrayfun(erk_constructor, ms, 'UniformOutput', false);
    methods        = [methods_nprk, methods_mrnprk, methods_erk];

    cb             = .2 * [1 1 1]; % black color for erk legend
    dc             = defaultColors();

    markers_nprk   = {ls('square--', dc(5,:), 5)};
    markers_mrnprk = {ls('<-', dc(1,:)), ls('>-', dc(2,:)), ls('*-', dc(3,:), 8), ls('diamond-', dc(4,:)), ls('o-', dc(5,:), 6), ls('^-', dc(6,:)), ls('hexagram-', dc(7,:))};
    markers_erk    = {ls('<:', cb),      ls('>:', cb),      ls('*:', cb, 8),      ls('diamond:', cb),      ls('o:', cb, 6),      ls('^:',cb),       ls('hexagram:', cb)};
    
    nm      = length(ms);
    markers = [markers_nprk, markers_mrnprk(1:nm), markers_erk(1:nm)];

    legend_order = [(2 : nm+1), 1, (nm + 2) : (2 * nm + 1)];

end

function s = ls(line_style, color, marker_size, line_width) % line style
    if(nargin < 3 || isempty(marker_size))
        marker_size = 4;
    end
    if(nargin < 4 || isempty(line_width))
        line_width = 4;
    end
    if(nargin < 2)
        color = [];
    end
    s = {line_style, 'LineWidth', line_width, 'MarkerSize', marker_size};
    if(~isempty(color))
        s = [s, {'Color', color}];
    end
end

function c = defaultColors()
    % result of get(gca, 'ColorOrder')
    c = [
             0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840
    ];

end

function f = makeConvergencePlot(tspan, Nts, error, markers, leg_str, plot_reorder)

    if(nargin < 6)
        plot_reorder = [];
    end

    f = figure();
    common_styles = {};
    
    for i = 1 : size(error,2)
        h = loglog(Nts, error(:, i), markers{i}{:}, common_styles{:}, 'MarkerIndices', 1 : 2 : length(Nts)); hold on;
        set(h, 'MarkerFaceColor', get(h,'Color'));
    end
    hold off;
    
    if(~isempty(leg_str))
        legend(leg_str, 'Location','northeastoutside'); legend box off;
    end
    if(~isempty(plot_reorder))
        chH = get(gca,'Children');
        set(gca, 'Children', [chH(plot_reorder)])
    end

    ylim([1e-12, 1e0]);
    xlim([Nts(1), Nts(end)]);
    ylabel(sprintf('Error at t=%.1f', tspan(end)), 'Interpreter', 'latex')
    xlabel('num steps $N_T$', 'Interpreter', 'latex')
    set(gca, 'FontSize', 14)
    set(gca,'TickLabelInterpreter','latex');
    grid on;
    grid minor;
    grid minor;
    
end

function f = makeEfficiencyPlot(tspan, time, error, leg_str, markers, plot_reorder)

    if(nargin < 6)
        plot_reorder = [];
    end

    f = figure();
    common_styles = {};
    
    for i = 1 : size(error,2)
        h = loglog(time(:,i), error(:, i), markers{i}{:}, common_styles{:}, 'MarkerIndices', 1 : 2 : size(time,1)); hold on;
        set(h, 'MarkerFaceColor', get(h,'Color'));
    end
    hold off;
    
    if(~isempty(leg_str))
        legend(leg_str, 'Location','northeastoutside'); legend box off;
    end
    if(~isempty(plot_reorder))
        chH = get(gca,'Children');
        set(gca, 'Children', [chH(plot_order)])
    end
    
    ylim([1e-11, 1e0]);
    xlim([5e-2, 1e2]);
    ylabel(sprintf('Error at t=%.1f', tspan(end)), 'Interpreter', 'latex')
    xlabel('runtime (s)', 'Interpreter', 'latex')
    set(gca, 'FontSize', 14);
    set(gca,'TickLabelInterpreter','latex');
    grid on;
    grid minor;
    grid minor;
    
end

function [f] = makeLegendPlot(markers, leg_str)

f = figure();
for i = 1 : length(markers)
    h = loglog(0, 0, markers{i}{:}); hold on;
    set(h, 'MarkerFaceColor', get(h,'Color')); 
end

legend(leg_str, 'Location', 'northeast', 'Orientation', 'horizontal', 'NumColumns', 4);
legend box off;

end

function [f] = makeSolutionPlot(prb)

[ts_plot, ys_plot] = prb.referenceSnapshots();

f = figure();
surf(prb.xs, ts_plot, ys_plot);
shading interp;
view([10 75]);
xlabel('x');
ylabel('t');
zlabel('u(x,t)');

end

function [integrators, line_styles, error, time] = reorderResults(integrators, line_styles, error, time, reorder)
    integrators = integrators(reorder);
    line_styles = line_styles(reorder);
    error       = error(:, reorder);
    time        = time(:, reorder);
end
