function dx = pa_model_mono(t, x, p, u_of_t)
% PA_MODEL_MONO  MONOSTABLE comparator (direct CYP11B2-turnover hypothesis).
% CLEAN VERSION: no salt-bridge feedback, and the contaminating fast-output
% gating phi(s) is REMOVED from the plastic part so this model cannot smuggle
% in fast-loop bistability. The only question this model answers is whether a
% NO-POSITIVE-FEEDBACK capacity dynamic can durably remit.
%
%   States [r; a; z_p; s]. Total capacity z = z_fixed + z_p.
%
%   Plastic capacity z_p:
%     'reversible'   (M1): z_p relaxes toward fixed disease set-point zp_set;
%                    drug adds direct decay -gu*u*z_p. No feedback, no
%                    irreversibility -> must relapse on withdrawal.
%     'irreversible' (M2): apoptosis ratchet -- the level z_p is pushed to does
%                    not recover (one-way). Biologically weakly supported (the
%                    monkey ZG apoptosis is largely salt-mediated and occurs
%                    amid net ZG GROWTH), included only as the strongest
%                    logically-possible competitor.
%
%   IMPORTANT (fix vs. earlier version): aldosterone OUTPUT here is NOT gated
%   by phi(s) for the plastic part. Output is linear in capacity. This prevents
%   the fast salt-output loop from creating a spurious low-output state that
%   masks capacity relapse. We measure remission on CAPACITY z_p, not a.

    r   = x(1); a = x(2); z_p = x(3);
    sodium = 0; if isfield(p,'sodium'), sodium = p.sodium; end
    S_target = a^p.ns/(p.Ks^p.ns + a^p.ns) + sodium; S_target = min(max(S_target,0),0.999);
    if isfield(p,'s_dynamic') && p.s_dynamic && numel(x)>=4, s = x(4); else, s = S_target; end

    z = p.zfixed + z_p;
    u = u_of_t(t);

    % Output linear in capacity (NO phi(s) gating) -> no fast-loop bistability.
    auto_prod  = p.alpha * (p.zfixed + z_p);
    renin_prod = p.beta * r^p.pr/(p.Kr^p.pr + r^p.pr);

    dr = p.r0*(1-s) - p.kr*r;
    da = (1-u)*(renin_prod + auto_prod) - p.ka*a;

    zset = p.zp_set;
    switch p.mono_mode
        case 'reversible'
            dz_p = p.eps*( -p.kz*(z_p - zset) - p.gu*u*z_p );
        case 'irreversible'
            eff_set = min(zset, z_p);                  % one-way ratchet
            dz_p = p.eps*( -p.kz*(z_p - eff_set) - p.gu*u*z_p );
        otherwise
            error('p.mono_mode must be reversible or irreversible');
    end

    if isfield(p,'s_dynamic') && p.s_dynamic && numel(x)>=4
        ds = (S_target - s)/p.tau_s;
        dx = [dr; da; dz_p; ds];
    else
        dx = [dr; da; dz_p];
    end
end
