function dx = pa_model(t, x, p, u_of_t)
% PA_MODEL  RAAS + adrenal-autonomy model for PA remission.
%
%   States (4-state version with dynamic salt signal):
%     r   = renin / AngII drive               (fast, ~hours)        x(1)
%     a   = circulating aldosterone           (fast, ~hours)        x(2)
%     z_p = PLASTIC autonomous capacity       (slow, ~weeks-months) x(3)
%     s   = salt / volume / MR state          (intermediate, ~weeks)x(4)
%
%   Total autonomous secretory capacity driving aldosterone production is
%       z = z_fixed + z_p
%   with:
%     z_fixed = IRREVERSIBLE (mutation-locked) autonomy. Pure parameter.
%               Secretes at FULL rate, UNGATED by the salt signal (truly
%               renin- and salt-independent) -> cannot be switched off by
%               output suppression. The mechanistic non-remitter floor.
%     z_p     = PLASTIC autonomy, salt-MODULATED, decays when salt support
%               is removed. The reversible part. (>= 0)
%
%   STEP 2 CHANGE -- dynamic salt signal s:
%     Previously s was algebraic: s = S(a) instantaneously. Now s is a STATE
%     with its own relaxation time tau_s:
%         tau_s * s_dot = S(a) - s,      S(a) = a^ns/(Ks^ns + a^ns)
%     This lags volume/salt/MR state behind circulating aldosterone. The point:
%     a SHORT drug pulse no longer instantly collapses the sustaining signal,
%     so crossing the separatrix requires SUSTAINED suppression -> a genuine
%     minimum treatment DURATION. Timescale ordering intended:
%         fast (r,a, ~hours)  <<  tau_s (~weeks)  <<  slow z_p (~months).
%     Set p.s_dynamic=0 to recover the old algebraic s=S(a) (consistency check).
%
%   Baxdrostat u(t) in [0,1] throttles aldosterone PRODUCTION only. It never
%   appears in the z_p or s equations.
%
%   t      : time (days)
%   x      : state vector [r; a; z_p; s]   (s optional; see below)
%   p      : parameter struct (see pa_params.m)
%   u_of_t : function handle u = u_of_t(t) in [0,1]
%   dx     : column vector of derivatives (same length as x)

    r   = x(1);
    a   = x(2);
    z_p = x(3);

    % Salt signal target: aldosterone-driven volume retention PLUS an additive
    % dietary-sodium contribution (sodium raises volume/salt state independent
    % of aldosterone). p.sodium in ~[0,0.5]; 0 = the original aldosterone-only
    % signal. Clamped to [0,1).
    sodium = 0; if isfield(p,'sodium'), sodium = p.sodium; end
    S_target = a^p.ns / (p.Ks^p.ns + a^p.ns) + sodium;
    S_target = min(max(S_target, 0), 0.999);           % keep in [0,1)
    if isfield(p,'s_dynamic') && p.s_dynamic && numel(x) >= 4
        s = x(4);
    else
        s = S_target;
    end

    z = p.zfixed + z_p;                          % total autonomous capacity
    u = u_of_t(t);                               % baxdrostat throttle in [0,1]

    % ----- Autonomous aldosterone production ----------------------------
    if p.fast_on
        actfrac = p.fmin + (1 - p.fmin) * s^p.nf / (p.Kf^p.nf + s^p.nf);
    else
        actfrac = 1;
    end
    if isfield(p,'legacy_z') && p.legacy_z
        % LEGACY (original 3-state) production: ALL capacity gated by phi(s),
        % single z decaying toward zfloor. Reproduces the pre-Step-1 results.
        auto_prod = p.alpha * z * actfrac;
    else
        % NEW: irreversible floor ungated + plastic part gated.
        auto_prod = p.alpha * ( p.zfixed + z_p * actfrac );
    end

    % ----- Renin-driven aldosterone production --------------------------
    renin_prod = p.beta * r^p.pr / (p.Kr^p.pr + r^p.pr);

    % ----- Fast subsystem: r and a --------------------------------------
    dr = p.r0 * (1 - s) - p.kr * r;
    da = (1 - u) * (renin_prod + auto_prod) - p.ka * a;

    % ----- Slow subsystem: PLASTIC capacity z_p -------------------------
    if p.slow_on
        sustain = p.gz * s^p.nz / (p.Kz^p.nz + s^p.nz);
    else
        sustain = 0;
    end
    if isfield(p,'legacy_z') && p.legacy_z
        % LEGACY: single z decays toward zfloor (here z==zfixed+z_p with
        % zfixed=0, so z_p decays toward zfloor).
        dz_p = p.eps * ( -p.kz * (z_p - p.zfloor) + sustain );
    else
        dz_p = p.eps * ( -p.kz * z_p + sustain );
    end

    % ----- Intermediate subsystem: salt/volume state s ------------------
    if isfield(p,'s_dynamic') && p.s_dynamic && numel(x) >= 4
        ds = ( S_target - s ) / p.tau_s;
        dx = [dr; da; dz_p; ds];
    else
        dx = [dr; da; dz_p];
    end
end
