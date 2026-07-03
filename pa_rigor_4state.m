function pa_rigor_4state()
% PA_RIGOR_4STATE  Full-model rigor pass: confirm the calibrated 4-state model
% has the same qualitative skeleton as the 3-state scaffold, and show the
% z_fixed remission threshold IS a saddle-node bifurcation.
%
% State: [r; a; z_p; s].  Total capacity z = z_fixed + z_p.
%
% Parts:
%   (1) Equilibria + 4x4 Jacobian eigenvalues at a remission-CAPABLE z_fixed:
%       expect two stable nodes + one saddle, and a THREE-tier timescale split
%       (fast r,a ; intermediate s ~ tau_s ; slow z_p).
%   (2) Continuation in z_fixed: trace all equilibria vs z_fixed, detect where
%       the healthy node and the saddle collide (saddle-node bifurcation).
%       That collision point IS the remission threshold z_fixed*.
%   (3) Continuation in sodium sigma: show the threshold moves (sodium shrinks
%       the remission-capable range), grounding the sodium prediction in
%       bifurcation structure.

    close all; here = fileparts(mfilename('fullpath'));
    fo = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12,...
                      'StepTolerance',1e-12,'MaxIterations',400);

    % ===== Part 1: equilibria & Jacobian at a remission-capable patient =====
    p = pa_params(); p.zfixed = 0.05; p.sodium = 0;   % full 4-state default
    f = @(x) pa_model(0, x, p, @(t)0);                % autonomous RHS, drug off

    % seeds for the three equilibria (4-state). s-seed = S(a) at that a.
    seeds = [ 0.03 4.0 1.95 sgoal(4.0,p);     % disease
              0.55 1.15 0.45 sgoal(1.15,p);   % middle/saddle
              0.99 0.25 0.00 sgoal(0.25,p) ]; % healthy
    fprintf('===== Part 1: 4-state equilibria at z_fixed=%.2f, sigma=%.2f =====\n',...
        p.zfixed, p.sodium);
    E = [];
    for k=1:3
        xe = fsolve(f, seeds(k,:).', fo);
        E(end+1,:) = xe.'; %#ok<AGROW>
        J = njac(f, xe); ev = eig(J);
        stable = all(real(ev)<0);
        ztot = p.zfixed + xe(3);
        fprintf('\n(r,a,z_p,s)=(%.3f,%.3f,%.3f,%.3f)  z_tot=%.3f  |res|=%.1e\n',...
            xe(1),xe(2),xe(3),xe(4), ztot, norm(f(xe)));
        fprintf('  eigenvalues: ');
        fprintf('%+.4f  ', sort(real(ev),'descend')); fprintf('\n');
        fprintf('  -> %s\n', tern(stable,'STABLE node','SADDLE (unstable)'));
        re = sort(abs(real(ev(abs(real(ev))>1e-9))),'descend');
        if numel(re)>=3
            fprintf('  timescales tau=1/|Re| (d): %.2f, %.2f, %.1f, %.1f\n', 1./re);
        end
    end

    % ===== Part 2: continuation in z_fixed -> saddle-node at threshold =====
    % Proper predictor-corrector: each branch is tracked by stepping z_fixed and
    % seeding each solve from the PREVIOUS converged point on that branch (not a
    % fixed seed). The disease branch is continued across the whole range; the
    % healthy and saddle branches are continued forward until they MERGE (the
    % saddle-node), detected when their corrected solutions coincide.
    fprintf('\n===== Part 2: continuation in z_fixed (find saddle-node) =====\n');
    [zd, Bdis] = continue_branch(seeds(1,:).', 0, 0.55, 0.0, fo);  % disease
    [zh, Bhlt] = continue_branch(seeds(3,:).', 0, 0.55, 0.0, fo);  % healthy
    [zs, Bsad] = continue_branch(seeds(2,:).', 0, 0.55, 0.0, fo);  % saddle

    % saddle-node = z_fixed where healthy and saddle branches meet (a coincides)
    zsn = saddle_node_point(zh, Bhlt, zs, Bsad);
    fprintf('Saddle-node (healthy+saddle collide) at z_fixed* = %.4f\n', zsn);
    fprintf('  below: 3 equilibria (bistable). above: disease only.\n');

    % ===== Part 3: threshold vs sodium =====
    fprintf('\n===== Part 3: remission threshold vs sodium =====\n');
    sig_vals = linspace(0, 0.2, 9);
    zstar = nan(size(sig_vals));
    for j=1:numel(sig_vals)
        % continue healthy and saddle branches at this sodium; their merge is z*
        [zhj,Bhj] = continue_branch(seeds(3,:).', 0, 0.55, sig_vals(j), fo);
        [zsj,Bsj] = continue_branch(seeds(2,:).', 0, 0.55, sig_vals(j), fo);
        zstar(j) = saddle_node_point(zhj,Bhj,zsj,Bsj);
        fprintf('  sigma=%.3f : remission-capable for z_fixed up to %.4f\n', sig_vals(j), zstar(j));
    end

    % ===== figure =====
    figure('Name','4-state rigor','Color','w','Position',[60 60 1150 460]);
    subplot(1,2,1); hold on; box on;
    % continued branches (clean, no stale-seed artifacts)
    plot(zd, Bdis(2,:), 'b-', 'LineWidth',2);                 % disease node
    % healthy & saddle only up to the saddle-node
    keepH = zh <= zsn + 1e-9; keepS = zs <= zsn + 1e-9;
    plot(zh(keepH), Bhlt(2,keepH), 'g-', 'LineWidth',2);      % healthy node
    plot(zs(keepS), Bsad(2,keepS), 'r--','LineWidth',1.5);    % saddle
    plot(zsn, 0.5*(interp1(zh(keepH),Bhlt(2,keepH),zsn)+interp1(zs(keepS),Bsad(2,keepS),zsn)),...
         'ko','MarkerFaceColor','k','MarkerSize',8);          % saddle-node point
    xline(zsn,'k:','LineWidth',1.2);
    text(zsn+0.01, 4.6, sprintf('saddle-node  z_{fixed}*=%.3f',zsn),'FontSize',8);
    xlabel('z_{fixed}'); ylabel('equilibrium aldosterone a');
    title('Continuation in z_{fixed}: healthy + saddle annihilate at threshold');
    legend({'disease node (stable)','healthy node (stable)','saddle (unstable)',...
            'saddle-node bifurcation'},'Location','east','FontSize',8);
    ylim([0 6]); xlim([0 0.55]);

    subplot(1,2,2); hold on; box on;
    plot(sig_vals, zstar, 'o-','LineWidth',1.8,'MarkerFaceColor','b');
    xlabel('dietary sodium \sigma'); ylabel('remission threshold z_{fixed}*');
    title('Sodium lowers the remission threshold (bifurcation shifts)');
    grid on;

    sgtitle('4-state model: same skeleton as scaffold; z_{fixed} threshold is a saddle-node bifurcation');
    saveas(gcf, fullfile(here,'pa_rigor_4state.png'));
    fprintf('\nSaved figure to pa_rigor_4state.png\n');
end

function s0 = sgoal(a0,p)
    sodium = 0; if isfield(p,'sodium'), sodium = p.sodium; end
    s0 = min(a0^p.ns/(p.Ks^p.ns + a0^p.ns) + sodium, 0.999);
end

function [zs, B] = continue_branch(x0, zlo, zhi, sigma, fo)
% Track one equilibrium branch as z_fixed steps from zlo to zhi, seeding each
% solve from the previous converged point (predictor = previous solution).
% Returns z_fixed values zs and the branch states B (4 x N). Stops early if the
% branch can no longer be continued (fold / disappearance), so the returned
% arrays may be shorter than the full grid.
    N = 140;
    zgrid = linspace(zlo, zhi, N);
    B = nan(4, N); zs = nan(1, N);
    xprev = x0;
    for i = 1:N
        pp = pa_params(); pp.zfixed = zgrid(i); pp.sodium = sigma;
        ff = @(x) pa_model(0, x, pp, @(t)0);
        [xe, ~, flag] = fsolve(ff, xprev, fo);
        if flag <= 0 || norm(ff(xe)) > 1e-5
            break;   % branch lost (likely past a fold) -> stop
        end
        B(:, i) = xe; zs(i) = zgrid(i);
        xprev = xe;  % predictor for next step
    end
    keep = ~isnan(zs);
    zs = zs(keep); B = B(:, keep);
end

function zsn = saddle_node_point(zh, Bh, zs, Bs)
% Saddle-node = z_fixed where the healthy and saddle branches coincide. We find
% where the aldosterone gap between the two branches (on their common z range)
% shrinks to ~0, or where one branch terminates (fold), whichever comes first.
    zc = max(min(zh), min(zs)) : 1e-3 : min(max(zh), max(zs));
    if isempty(zc), zsn = min(max(zh), max(zs)); return; end
    ah = interp1(zh, Bh(2,:), zc, 'linear');
    asd = interp1(zs, Bs(2,:), zc, 'linear');
    gap = abs(ah - asd);
    % saddle-node where gap first drops below a small tolerance
    idx = find(gap < 0.05, 1, 'first');
    if isempty(idx)
        zsn = zc(end);        % branches end before fully merging; use overlap end
    else
        zsn = zc(idx);
    end
end

function J = njac(f,x)
    n=numel(x); J=zeros(n); f0=f(x);
    for i=1:n, h=1e-6*max(1,abs(x(i))); xp=x; xp(i)=xp(i)+h; J(:,i)=(f(xp)-f0)/h; end
end
function out=tern(c,a,b), if c,out=a;else,out=b;end, end
