function pa_bifurcation_sensitivity()
% PA_BIFURCATION_SENSITIVITY  Local parameter sensitivity of the saddle-node
% threshold z_fixed* (Gemini reviewer point 2). Fast version.
%
% For each key parameter: perturb +/-10% and +/-20% (one at a time), recompute
% z_fixed*(sigma=0), report (i) whether the saddle-node still EXISTS and
% (ii) how far z_fixed* MOVES.
%
% eps (slow timescale) and tau_s (salt lag) are CONTROLS: they set speed, not
% equilibrium structure, so z_fixed* should be ~invariant to them.

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',5.0);
    Tset = 3000;   % settle time (~30 slow time-constants); plenty for equilibrium

    base = pa_params();
    z0 = find_threshold(base, opts, Tset);
    fprintf('\nBaseline z_fixed* (sigma=0) = %.4f\n\n', z0);

    params = {'gz','kz','Kz','nz','Ks','ns','alpha','beta','tau_s','eps'};
    pct    = [-20 -10 10 20];

    fprintf('%-8s %8s', 'param', 'base');
    for q=pct, fprintf('  %+5d%%', q); end
    fprintf('   exists\n');
    fprintf('%s\n', repmat('-',1,64));

    maxmove=0; maxmove_p='';
    for i=1:numel(params)
        pname=params{i}; bval=base.(pname);
        fprintf('%-8s %8.3f', pname, bval);
        allexist=true;
        for q=pct
            p=base; p.(pname)=bval*(1+q/100);
            zc=find_threshold(p,opts,Tset);
            if isnan(zc)||zc<=0.021
                fprintf('    n/a '); allexist=false;
            else
                fprintf('   %.3f', zc);
                mv=abs(zc-z0)/z0*100;
                if mv>maxmove, maxmove=mv; maxmove_p=sprintf('%s@%+d%%',pname,q); end
            end
        end
        fprintf('   %s\n', tern(allexist,'yes','CHECK'));
    end
    fprintf('%s\n', repmat('-',1,64));
    fprintf('Largest relative shift in z_fixed* over all +/-20%% perturbations: %.1f%% (%s)\n', maxmove, maxmove_p);
    fprintf('(eps, tau_s are structural controls: expected ~0%% shift.)\n');
end

function zc = find_threshold(p, opts, Tset)
    zc=NaN; lastgood=NaN;
    for zf=0.02:0.01:1.3
        pp=p; pp.zfixed=zf;
        [~,X]=ode15s(@(t,x) pa_model(t,x,pp,@(tt)0),[0 Tset],[0.98;0.34;0.0;0.02],opts);
        if X(end,2)<1.0, lastgood=zf;
        else, if ~isnan(lastgood), zc=0.5*(lastgood+zf); return; end
        end
    end
    zc=lastgood;
end

function o=tern(c,a,b), if c, o=a; else, o=b; end, end
