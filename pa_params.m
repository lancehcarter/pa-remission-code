function p = pa_params()
% PA_PARAMS  Default parameters for the PA remission scaffold model.
%
% Timescale note: time is in DAYS. Fast variables (r,a) have rates ~O(1-10)/day
% so they equilibrate in well under a day (consistent with hours). The slow
% variable z is driven by p.eps << 1, so it evolves over weeks-months.
%
% These values are a STARTING POINT chosen to put the system in a regime where
% bistability is *possible* but not hand-forced. They are not fitted to data.
% The whole point of the scaffold is to vary these and see what the dynamics do.

    % ---- Fast aldosterone (a) ----
    % NB: rates are in nondimensional production units, not literal hormone
    % half-lives. What matters dynamically is that a,r are fast relative to z
    % (ensured by p.eps), not the absolute value of ka. ka was rescaled from an
    % initial (too-large) literal half-life value that over-damped the system
    % and destroyed the high-aldosterone disease state.
    p.ka    = 1.0;     % aldosterone clearance rate (fast relative to z)
    p.alpha = 2.0;     % autonomous secretion gain (per unit capacity)
    p.beta  = 0.5;     % renin-driven secretion gain
    p.Kr    = 1.0;     % half-sat of renin drive in aldosterone production
    p.pr    = 2;       % Hill coeff, renin -> aldosterone

    % ---- Fast renin (r) ----
    p.r0    = 4.0;     % renin baseline drive (1/day)
    p.kr    = 4.0;     % renin relaxation rate (1/day)

    % ---- Salt / volume signal s(a) ----
    p.Ks    = 1.2;     % half-sat: aldosterone level at which volume signal is half-max
    p.ns    = 3;       % steepness of volume response to aldosterone

    % ---- STEP 2: dynamic salt signal (lag) ----
    % s becomes a state with relaxation time tau_s: tau_s*s_dot = S(a) - s.
    % Intermediate timescale: slower than fast (r,a ~hours), faster than slow
    % z_p (~100 d). ~2 weeks for volume/sodium/MR normalization. Set
    % s_dynamic=0 to recover the algebraic s=S(a).
    p.s_dynamic = 1;   % 1 = s is a dynamic state with lag; 0 = algebraic
    p.tau_s = 14;      % salt/volume relaxation time (days)

    % ---- STEP 3: dietary sodium (additive offset to salt signal) ----
    % Raises the salt/volume target independent of aldosterone. Higher sodium
    % -> stronger sustain of z_p -> harder to remit. This is the handle for the
    % testable sodium-dependence-of-remission prediction. Range ~[0, 0.4].
    p.sodium = 0;      % dietary sodium contribution to salt signal

    % ---- FAST self-reinforcement (output bistability candidate) ----
    p.fast_on = 1;     % 1 = fast salt/volume -> autonomous-output loop ON
    p.fmin  = 0.05;    % minimum active fraction (floor when volume signal low)
    p.Kf    = 0.45;    % half-sat of active-fraction sigmoid in s
    p.nf    = 6;       % steepness (high => switch-like fast output)

    % ---- SLOW capacity: plastic part z_p (salt-bridge bistability) ----
    p.slow_on = 1;     % 1 = salt bridge ON (z_p sustained by salt/volume state)
    p.eps   = 0.01;    % slow timescale separation (z_p evolves ~100x slower)
    p.kz    = 1.0;     % plastic-capacity decay rate (scaled by eps)
    p.gz    = 2.0;     % strength of salt/volume sustaining z_p
    p.Kz    = 0.5;     % half-sat of sustaining sigmoid in s
    p.nz    = 6;       % steepness of sustaining sigmoid
    p.legacy_z = 0;    % 0 = new z_p dynamics; 1 = original single-z (see pa_params_legacy)
    p.zfloor = 0.05;   % legacy capacity floor (only used when legacy_z=1)

    % ---- IRREVERSIBLE capacity z_fixed (per-patient; the key new knob) ----
    % Total capacity z = z_fixed + z_p. The drug and salt bridge act ONLY on
    % z_p; z_fixed is mutation-locked and untouchable. If z_fixed sits above
    % the saddle, the patient CANNOT remit by output suppression alone.
    % Default is low (below saddle) so the baseline patient is remission-capable;
    % this parameter is swept to map who can remit.
    p.zfixed = 0.05;   % irreversible autonomous capacity floor (per-patient)

    % ---- Disease vs healthy reference capacities (for initial conditions) ----
    p.z_disease = 2.0; % high TOTAL autonomous capacity (PA state)
    p.z_health  = p.zfixed;  % low capacity (remitted): plastic part gone, only floor
end
