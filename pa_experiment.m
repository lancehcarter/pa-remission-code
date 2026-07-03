function pa_experiment()
% PA_EXPERIMENT  Corrected scaffold experiment (pure-coupling version).
%
% For each mechanism scenario we ask THREE things in order:
%   (1) EXISTENCE: drug-free, does a DISEASE equilibrium exist at all?
%       (settle from a high start and a low start; do they differ?)
%   (2) If a disease state exists, apply a FINITE baxdrostat pulse starting
%       FROM that disease equilibrium, then withdraw.
%   (3) CLASSIFY: remission only if we started at a genuine disease state and
%       ended at the healthy (low) equilibrium after withdrawal.
%
% This avoids the earlier artifact where z simply decayed and every scenario
% trivially "remitted". Now "neither" should correctly report NO DISEASE STATE.

    close all;
    here = fileparts(mfilename('fullpath'));

    % ---- Drug pulse (days) ----
    t_start = 50;
    t_stop  = 50 + 16*30;     % ~16 months
    u_max   = 0.95;
    u_of_t  = @(t) u_max * (t >= t_start && t < t_stop);
    t_end   = t_stop + 720;   % 2 yr follow-up after withdrawal

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1.0);
    u0   = @(t) 0;            % drug-free

    scenarios = { ...
        'A: both layers', 1, 1; ...
        'B: fast only',   1, 0; ...
        'C: slow only',   0, 1; ...
        'D: neither',     0, 0  };

    figure('Name','PA remission scaffold','Color','w','Position',[100 80 1100 850]);

    fprintf('\n%-16s | disease eq (a,z) | healthy eq (a,z) | bistable? | pulse outcome\n', 'scenario');
    fprintf('%s\n', repmat('-',1,92));

    for i = 1:size(scenarios,1)
        name = scenarios{i,1};
        p = pa_params_legacy(); p.fast_on = scenarios{i,2}; p.slow_on = scenarios{i,3};

        % (1) EXISTENCE: settle drug-free from HIGH and LOW starts.
        xh = settle(@(t,x) pa_model(t,x,p,u0), [0.05; 4.0; p.z_disease], opts);
        xl = settle(@(t,x) pa_model(t,x,p,u0), [1.5;  0.2; p.zfloor],    opts);
        a_hi = xh(2); z_hi = xh(3);
        a_lo = xl(2); z_lo = xl(3);
        bistable = abs(a_hi - a_lo) > 0.25*max(a_hi,1e-6);   % distinct states?

        if ~bistable
            fprintf('%-16s | (only one state: a=%.3f z=%.3f)            | %-9s | %s\n', ...
                name, a_hi, z_hi, 'NO', 'N/A (no disease state to remit from)');
            outcome = 'no disease state';
            % still simulate from high start for the plot
            x0 = xh;
        else
            % (2) PULSE from the DISEASE equilibrium.
            x0 = xh;
        end

        [t, X] = ode15s(@(t,x) pa_model(t,x,p,u_of_t), [0 t_end], x0, opts);
        a = X(:,2); z = X(:,3);
        uvec = arrayfun(u_of_t, t);

        if bistable
            a_final = mean(a(t > t_end - 180));
            % remission = ended near the healthy low state, well below disease
            remitted = (a_final < 0.5*a_hi) && (abs(a_final - a_lo) < 0.5*max(a_lo,0.2));
            outcome = ternary(remitted, 'DURABLE REMISSION', 'relapse (returned to disease)');
            fprintf('%-16s | a=%.3f z=%.3f      | a=%.3f z=%.3f     | %-9s | %s\n', ...
                name, a_hi, z_hi, a_lo, z_lo, 'YES', outcome);
        end

        % ---- Plot ----
        subplot(size(scenarios,1),1,i);
        yyaxis left
        plot(t, a, '-', 'LineWidth', 1.6); hold on;
        plot(t, z, '--', 'LineWidth', 1.1);
        ylabel('a (—), z (--)');
        ylim([0 max(4.2, p.z_disease*1.1)]);
        yyaxis right
        area(t, uvec, 'FaceAlpha', 0.10, 'EdgeColor','none');
        ylabel('drug u'); ylim([0 1.2]);
        xline(t_start,':'); xline(t_stop,':');
        title(sprintf('%s   ->   %s', name, outcome));
        xlabel('time (days)');
    end

    sgtitle('Finite baxdrostat pulse from the disease equilibrium (pure-coupling model)');
    saveas(gcf, fullfile(here,'pa_experiment.png'));
    fprintf('\nSaved figure to pa_experiment.png\n');
end

% ---- helpers ----
function xend = settle(odefun, x0, opts)
    [~, X] = ode15s(odefun, [0 8000], x0, opts);
    xend = X(end,:).';
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
