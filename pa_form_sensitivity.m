function pa_form_sensitivity()
% PA_FORM_SENSITIVITY  Is the slow bistability an artifact of the Hill shape?
%
% The reduced slow dynamics are  zdot = eps * F(z),  with
%     F(z) = -kz*(z - zfloor) + gz * B( s(A(z)) )
% where:
%   A(z) = fast-subsystem aldosterone equilibrium given frozen capacity z
%   s(a) = a^ns/(Ks^ns+a^ns)             (salt/volume signal)
%   B(.) = the SALT-BRIDGE shape under test (input in [0,1], output in [0,1])
%
% Bistability  <=>  F(z) has 3 zeros (two stable, one unstable),
% equivalently the sustain curve gz*B(s(A(z))) crosses the decay line
% kz*(z-zfloor) three times.
%
% We swap B for several increasing, saturating shapes and ask whether
% the 3-crossing structure survives, and over how wide a gz range.

    close all;
    here = fileparts(mfilename('fullpath'));
    p = pa_params_legacy();
    p.slow_on = 1;            % bridge present (its SHAPE is what we vary)
    p.fast_on = 1;            % keep fast loop on (we showed it's ~inert; check holds)

    % ---- candidate bridge shapes B(s), all increasing, B(0)=0-ish, B(1)~1 ----
    nz = p.nz; Kz = p.Kz;
    shapes = struct( ...
        'name', {'Hill (n=6)', 'Hill (n=2)', 'Michaelis (n=1)', ...
                 'saturating exp', 'piecewise-linear', 'tanh-sigmoid'}, ...
        'B', { @(s) s.^6 ./(Kz^6 + s.^6), ...
               @(s) s.^2 ./(Kz^2 + s.^2), ...
               @(s) s    ./(Kz   + s), ...
               @(s) 1 - exp(-s/Kz), ...
               @(s) min(max((s-0.2)/0.4,0),1), ...
               @(s) 0.5*(1+tanh((s-Kz)/0.12)) } );

    % ---- slow manifold A(z): fast-subsystem aldosterone equilibrium ----
    Afun = @(z) fast_equilibrium_a(z, p);

    zgrid = linspace(0, 2.5, 600);
    Avals = arrayfun(Afun, zgrid);
    svals = Avals.^p.ns ./ (p.Ks^p.ns + Avals.^p.ns);
    decay = p.kz*(zgrid - p.zfloor);

    figure('Name','Salt-bridge form sensitivity','Color','w','Position',[80 80 1150 720]);

    % Panel 1: the crossing geometry for each shape (at default gz)
    subplot(1,2,1); hold on; box on;
    plot(zgrid, decay, 'k-', 'LineWidth', 2, 'DisplayName','decay line');
    nshapes = numel(shapes);
    cols = lines(nshapes);
    fprintf('\nShape                | #equilibria (default gz=%.2f) | bistable?\n', p.gz);
    fprintf('%s\n', repmat('-',1,68));
    for i = 1:nshapes
        sustain = p.gz * shapes(i).B(svals);
        plot(zgrid, sustain, '-', 'Color', cols(i,:), 'LineWidth', 1.4, ...
             'DisplayName', shapes(i).name);
        F = sustain - decay;
        nz_cross = sum(F(1:end-1).*F(2:end) < 0);
        fprintf('%-20s | %-30d | %s\n', shapes(i).name, nz_cross, ...
            ternary(nz_cross>=3,'YES','no'));
    end
    xlabel('capacity z'); ylabel('sustain  gz·B(s(A(z)))   vs   decay');
    title(sprintf('Crossing geometry at gz=%.2f', p.gz));
    legend('Location','northwest'); 

    % Panel 2: for each shape, the RANGE of gz giving bistability
    subplot(1,2,2); hold on; box on;
    gzrange = linspace(0.2, 6, 200);
    for i = 1:nshapes
        bist = false(size(gzrange));
        for k = 1:numel(gzrange)
            sustain = gzrange(k) * shapes(i).B(svals);
            F = sustain - decay;
            bist(k) = sum(F(1:end-1).*F(2:end) < 0) >= 3;
        end
        % plot a band at height i where bistable
        yy = i*ones(size(gzrange)); yy(~bist) = NaN;
        plot(gzrange, yy, '-', 'Color', cols(i,:), 'LineWidth', 7);
    end
    set(gca,'YTick',1:nshapes,'YTickLabel',{shapes.name});
    xlabel('salt-bridge strength  g_z'); ylabel('bridge shape');
    title('Range of g_z for which bistability exists');
    ylim([0.5 nshapes+0.5]); xlim([gzrange(1) gzrange(end)]);

    sgtitle('Does slow bistability survive changes in the salt-bridge functional form?');
    saveas(gcf, fullfile(here,'pa_form_sensitivity.png'));
    fprintf('\nSaved figure to pa_form_sensitivity.png\n');
end

% ---- fast subsystem aldosterone equilibrium at frozen z (drug off) ----
function aeq = fast_equilibrium_a(z, p)
    % Solve k_a*a = renin_prod(r_eq(a)) + alpha*z*phi(s(a)) for a, on the
    % UPPER branch (disease side) when multiple roots exist; else the unique root.
    avals = linspace(0, 10, 4000);
    g = zeros(size(avals));
    for i = 1:numel(avals)
        a = avals(i);
        s = a^p.ns/(p.Ks^p.ns + a^p.ns);
        r = p.r0*(1-s)/p.kr;
        if p.fast_on
            phi = p.fmin + (1-p.fmin)*s^p.nf/(p.Kf^p.nf + s^p.nf);
        else
            phi = 1;
        end
        rp = p.beta*r^p.pr/(p.Kr^p.pr + r^p.pr);
        g(i) = rp + p.alpha*z*phi - p.ka*a;     % = 0 at fast equilibria
    end
    sc = find(g(1:end-1).*g(2:end) < 0);
    if isempty(sc)
        [~,im] = min(abs(g)); aeq = avals(im);
    else
        % take the largest root (upper/disease branch) for the manifold we ride
        aeq = avals(sc(end));
    end
end

function out = ternary(cond,a,b), if cond, out=a; else, out=b; end, end
