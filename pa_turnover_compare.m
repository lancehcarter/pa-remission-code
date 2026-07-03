function pa_turnover_compare()
% PA_TURNOVER_COMPARE  Fair-comparator test (addresses reviewer point #1).
%
% Compares the BISTABLE salt-bridge model against a two-compartment DIRECT-
% TURNOVER model that has patient-specific resistant (drug-insensitive) and
% sensitive (drug-driven, irreversible) capacity compartments but NO feedback
% and NO bistability.
%
% Both models are given the SAME patient axis: the absolute amount of
% irreversible / resistant autonomous capacity (z_fixed for the bistable model,
% z_R for the turnover model). For each value we apply a long full-dose course,
% withdraw, settle to the true drug-free asymptote, and record post-withdrawal
% aldosterone.
%
% THE QUESTION: does durable remission appear as a SHARP threshold (a
% discontinuity / bifurcation) or a SMOOTH gradient in the patient parameter?
%   - Bistable model: expect a discontinuity at the saddle-node z_fixed* (the
%     healthy equilibrium ceases to exist).
%   - Turnover model: expect a SMOOTH, continuous boundary (post-withdrawal
%     aldosterone is a continuous function of z_R; no equilibrium vanishes).
% A discontinuous vs continuous boundary is a genuine qualitative discriminator
% and a testable difference, not an artifact of a rigged competitor.

    close all; here = fileparts(mfilename('fullpath'));
    opts = odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',2.0);
    t_pre=50; t_treat=24*30; t_post=200; u_max=0.95;   % long course (2 yr)

    zgrid = linspace(0, 0.8, 81);     % common patient axis: resistant capacity

    % ---------- BISTABLE model: sweep z_fixed ----------
    aB = nan(size(zgrid)); aB_dis = nan(size(zgrid));
    for i=1:numel(zgrid)
        p = pa_params(); p.zfixed = zgrid(i);
        s0 = min(4^p.ns/(p.Ks^p.ns+4^p.ns),0.999);
        xeq = settle(@pa_model, [0.05;4;max(p.z_disease-p.zfixed,0.05);s0], p, opts);
        aB_dis(i) = xeq(2);
        uft = @(t) u_max*(t>=t_pre && t<t_pre+t_treat);
        [t,X] = ode15s(@(t,x)pa_model(t,x,p,uft), [0 t_pre+t_treat+t_post], xeq, opts);
        xf = settle(@pa_model, X(end,:).', p, opts);
        aB(i) = xf(2);
    end

    % ---------- TURNOVER model: sweep z_R ----------
    aT = nan(size(zgrid)); aT_dis = nan(size(zgrid));
    z_tot = 2.0;                       % same total disease capacity reference
    for i=1:numel(zgrid)
        p = pa_params(); p = turnover_params(p);
        p.zR = zgrid(i);
        zS0 = max(z_tot - p.zR, 0);     % sensitive = remainder of disease capacity
        s0 = min(4^p.ns/(p.Ks^p.ns+4^p.ns),0.999);
        xeq = settle(@pa_model_turnover, [0.05;4;zS0;s0], p, opts);
        aT_dis(i) = xeq(2);
        uft = @(t) u_max*(t>=t_pre && t<t_pre+t_treat);
        [t,X] = ode15s(@(t,x)pa_model_turnover(t,x,p,uft), [0 t_pre+t_treat+t_post], xeq, opts);
        xf = settle(@pa_model_turnover, X(end,:).', p, opts);
        aT(i) = xf(2);
    end

    % ---------- sharpness metric: max |d a_final / d z| ----------
    gB = max(abs(diff(aB)./diff(zgrid)));
    gT = max(abs(diff(aT)./diff(zgrid)));
    % remission thresholds (where post-withdrawal a first exceeds half disease)
    thrB = first_cross(zgrid, aB, 0.5*max(aB_dis));
    thrT = first_cross(zgrid, aT, 0.5*max(aT_dis));

    fprintf('\n--- Fair comparator: bistable vs two-compartment turnover ---\n');
    fprintf('BISTABLE : remission threshold z_fixed* = %.3f ; max slope |da/dz| = %.1f\n', thrB, gB);
    fprintf('TURNOVER : remission threshold z_R*     = %.3f ; max slope |da/dz| = %.1f\n', thrT, gT);
    fprintf('Slope ratio (bistable/turnover) = %.1f  (>>1 means bistable boundary is much sharper)\n', gB/gT);

    % ---------- plot ----------
    figure('Name','Fair comparator','Color','w','Position',[60 60 1150 640]);

    subplot(2,2,[1 3]);
    plot(zgrid, aB, 'LineWidth', 2.0); hold on;
    plot(zgrid, aT, 'LineWidth', 2.0);
    if ~isnan(thrB), xline(thrB,'--','Color',[0 0.45 0.74]); end
    if ~isnan(thrT), xline(thrT,'--','Color',[0.85 0.33 0.10]); end
    xlabel('irreversible / resistant capacity  (z_{fixed} or z_R)');
    ylabel('post-withdrawal aldosterone  a_{final}');
    legend({'bistable (salt bridge)','two-compartment turnover'},'Location','southeast');
    title('Post-withdrawal aldosterone vs patient-intrinsic resistant capacity');
    grid on; ylim([0 4.5]);

    subplot(2,2,2);
    % numerical derivative (boundary sharpness)
    zc = 0.5*(zgrid(1:end-1)+zgrid(2:end));
    plot(zc, abs(diff(aB)./diff(zgrid)),'LineWidth',1.6); hold on;
    plot(zc, abs(diff(aT)./diff(zgrid)),'LineWidth',1.6);
    xlabel('z_{fixed} or z_R'); ylabel('|d a_{final}/dz|');
    title('Boundary sharpness (derivative)');
    legend({'bistable','turnover'},'Location','northeast'); grid on;

    subplot(2,2,4);
    plot(zgrid, aB_dis,'LineWidth',1.4); hold on;
    plot(zgrid, aT_dis,'LineWidth',1.4);
    xlabel('z_{fixed} or z_R'); ylabel('untreated disease a');
    title('Untreated disease aldosterone (context)');
    legend({'bistable','turnover'},'Location','best'); grid on; ylim([0 4.5]);

    sgtitle('Fair comparator: does remission appear as a sharp bifurcation or a smooth gradient?');
    saveas(gcf, fullfile(here,'pa_turnover_compare.png'));
    fprintf('\nSaved figure to pa_turnover_compare.png\n');
end

function p = turnover_params(p)
    p.slow_on = 0; p.fast_on = 0; p.legacy_z = 0;
    p.gu = 1.5;                 % drug-driven turnover gain (matches mono comparator)
    p.zR = 0.05;
end
function xe = settle(mdl, x0, p, opts)
    [~,X] = ode15s(@(t,x) mdl(t,x,p,@(tt)0), [0 12000], x0, opts);
    xe = X(end,:).';
end
function z = first_cross(xg, yg, level)
    idx = find(yg >= level, 1, 'first');
    if isempty(idx) || idx==1, z = NaN; else, z = xg(idx); end
end
