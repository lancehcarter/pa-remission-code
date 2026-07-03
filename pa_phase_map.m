function pa_phase_map()
% PA_PHASE_MAP  Step 3: the remission "prediction engine" phase diagram.
%
% Main panel: remission outcome over (treatment duration, z_fixed). This shows
% who remits and how long it takes, jointly. Expected structure:
%   - an UPPER bound on z_fixed (patient-intrinsic threshold; duration can't help)
%   - a LEFT/minimum-duration boundary that RISES as z_fixed approaches threshold
%     (near-threshold patients need longer treatment) -- the tau_s x z_fixed
%     coupling predicted in Step 2.
%
% Extra panels: how the remission region shifts with dietary SODIUM (the
% testable prediction) and with drug efficacy u_max.

    close all;
    here = fileparts(mfilename('fullpath'));

    t_pre = 50; t_post = 1000;
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',2.0);

    durations = round(logspace(log10(7), log10(720), 26));   % 1 wk .. 2 yr
    zfix_vals = linspace(0, 0.6, 26);

    % ---- main panel: baseline sodium=0, u_max=0.95 ----
    R_main = remission_grid(durations, zfix_vals, 0.0, 0.95, t_pre, t_post, opts);

    % ---- sodium comparison: low vs high sodium ----
    R_naLo = remission_grid(durations, zfix_vals, 0.00, 0.95, t_pre, t_post, opts);
    R_naHi = remission_grid(durations, zfix_vals, 0.15, 0.95, t_pre, t_post, opts);

    % ---- u_max comparison ----
    R_uHi  = remission_grid(durations, zfix_vals, 0.0, 0.95, t_pre, t_post, opts);
    R_uLo  = remission_grid(durations, zfix_vals, 0.0, 0.80, t_pre, t_post, opts);

    % ================= plotting =================
    figure('Name','PA remission phase map','Color','w','Position',[60 60 1250 780]);

    subplot(2,2,1);
    draw_map(durations, zfix_vals, R_main);
    title('Remission map  (sodium=0, u_{max}=0.95)');
    xlabel('treatment duration (days)'); ylabel('z_{fixed} (irreversible capacity)');

    subplot(2,2,2);
    % overlay sodium boundaries
    draw_map(durations, zfix_vals, R_naLo);
    hold on; contour(durations, zfix_vals, R_naHi, [0.5 0.5], 'm-', 'LineWidth', 2);
    title('Sodium effect: green=remit (low Na); magenta line = remission boundary (high Na)');
    xlabel('treatment duration (days)'); ylabel('z_{fixed}');
    text(0.05,0.95,'higher sodium shrinks remission region','Units','normalized',...
        'VerticalAlignment','top','FontSize',8,'BackgroundColor','w');

    subplot(2,2,3);
    draw_map(durations, zfix_vals, R_uHi);
    hold on; contour(durations, zfix_vals, R_uLo, [0.5 0.5], 'b-', 'LineWidth', 2);
    title('Drug efficacy: green=remit (u_{max}=0.95); blue line = boundary (u_{max}=0.80)');
    xlabel('treatment duration (days)'); ylabel('z_{fixed}');

    % summary panel: remitting AREA fraction vs sodium
    subplot(2,2,4); hold on; box on;
    na_vals = linspace(0, 0.25, 9);
    frac = zeros(size(na_vals));
    for i=1:numel(na_vals)
        R = remission_grid(durations, zfix_vals, na_vals(i), 0.95, t_pre, t_post, opts);
        frac(i) = mean(R(:));
    end
    plot(na_vals, 100*frac, 'o-', 'LineWidth',1.8, 'MarkerFaceColor','b');
    xlabel('dietary sodium offset'); ylabel('% of (duration,z_{fixed}) grid that remits');
    title('Testable prediction: remission likelihood falls with sodium');
    grid on;

    sgtitle('Step 3: remission phase map -- who remits, how long it takes, and sodium dependence');
    saveas(gcf, fullfile(here,'pa_phase_map.png'));

    % ---- text summary ----
    fprintf('Main map (sodium=0, u_max=0.95): %.0f%% of grid remits\n', 100*mean(R_main(:)));
    fprintf('Sodium 0.00 -> 0.15: remit fraction %.0f%% -> %.0f%%\n', ...
        100*mean(R_naLo(:)), 100*mean(R_naHi(:)));
    fprintf('u_max 0.95 -> 0.80: remit fraction %.0f%% -> %.0f%%\n', ...
        100*mean(R_uHi(:)), 100*mean(R_uLo(:)));
    fprintf('\nSaved figure to pa_phase_map.png\n');
end

% ===== run a full (duration x z_fixed) grid, return logical remission matrix =====
function R = remission_grid(durations, zfix_vals, sodium, umax, t_pre, t_post, opts)
    R = zeros(numel(zfix_vals), numel(durations));
    for iz = 1:numel(zfix_vals)
        p = pa_params();
        p.zfixed = zfix_vals(iz);
        p.sodium = sodium;
        % disease equilibrium (drug off) for this patient
        u0 = @(t) 0;
        x0hi = init_state(p, 0.05, 4.0, max(p.z_disease - p.zfixed, 0.1));
        xeq = settle(@(t,x) pa_model(t,x,p,u0), x0hi, opts);
        a_dis = xeq(2);
        % healthy reference
        xlo = settle(@(t,x) pa_model(t,x,p,u0), init_state(p,1.5,0.2,0.0), opts);
        a_hlt = xlo(2);
        disease_exists = abs(a_dis - a_hlt) > 0.25*max(a_dis,1e-6);
        for id = 1:numel(durations)
            if ~disease_exists
                R(iz,id) = 0;       % no disease state to remit from / monostable disease
                continue;
            end
            d = durations(id);
            uft = @(t) umax*(t >= t_pre && t < t_pre + d);
            [t,X] = ode15s(@(t,x) pa_model(t,x,p,uft), [0 t_pre+d+t_post], xeq, opts);
            af = mean(X(t > t(end)-150, 2));
            R(iz,id) = (af < 0.5*a_dis) && (abs(af - a_hlt) < 0.6*max(a_hlt,0.3));
        end
    end
end

function draw_map(durations, zfix_vals, R)
    imagesc(durations, zfix_vals, R); set(gca,'YDir','normal');
    colormap([0.95 0.85 0.85; 0.3 0.75 0.3]);   % red-ish = no, green = remit
    set(gca,'XScale','log'); xlim([durations(1) durations(end)]);
    caxis([0 1]);
end

function x0 = init_state(p, r0, a0, zp0)
    if isfield(p,'s_dynamic') && p.s_dynamic
        sodium = 0; if isfield(p,'sodium'), sodium = p.sodium; end
        s0 = min(a0^p.ns/(p.Ks^p.ns + a0^p.ns) + sodium, 0.999);
        x0 = [r0; a0; zp0; s0];
    else
        x0 = [r0; a0; zp0];
    end
end

function xend = settle(odefun, x0, opts)
    [~,X] = ode15s(odefun, [0 9000], x0, opts);
    xend = X(end,:).';
end
