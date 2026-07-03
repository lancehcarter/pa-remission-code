function p = pa_params_mono(mode)
% PA_PARAMS_MONO  Parameters for the monostable comparator (direct-turnover).
% Inherits the shared parameters from pa_params and adds the direct-turnover
% knobs. mode = 'reversible' (M1) or 'irreversible' (M2).
%
% Design intent: give the monostable hypothesis its BEST shot. The disease
% set-point zp_set is chosen so the untreated disease aldosterone matches the
% bistable model's disease state (~4), and the drug-turnover gain gu is set so
% that a long treatment meaningfully suppresses z_p. We then ASK whether this
% model can reproduce rarity + durability + minimum-duration jointly.

    if nargin < 1, mode = 'reversible'; end
    p = pa_params();           % shared fast/salt/zfixed parameters

    p.mono_mode = mode;
    p.slow_on = 0;             % no salt bridge (irrelevant in mono model anyway)
    p.legacy_z = 0;

    % disease maintenance set-point for plastic capacity (no-feedback baseline).
    % Matched so z_fixed + zp_set ~ z_disease (=2.0) at zfixed=0.05.
    p.zp_set = 1.95;

    % direct drug-driven turnover gain (CYP11B2 inhibition -> involution).
    % Chosen so that full-dose treatment over ~months substantially reduces z_p.
    p.gu = 1.5;
end
