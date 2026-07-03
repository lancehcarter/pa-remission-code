function p = pa_params_legacy()
% PA_PARAMS_LEGACY  Original (pre-calibration) parameter set.
%
% Returns the default parameters with ALL new modeling features turned OFF, so
% the original 3-state experiments (Experiments 1-5 in the development
% document) reproduce exactly:
%     legacy_z  = 1   -> single z decaying toward zfloor, fully salt-gated
%                        production (the original z dynamics)
%     s_dynamic = 0   -> algebraic salt signal s = S(a) (no lag)
%     zfixed    = 0   -> no separate irreversible floor (the floor is zfloor,
%                        carried inside the legacy z dynamics)
%     sodium    = 0   -> no dietary-sodium offset
%
% In legacy mode the state is 3-dimensional [r; a; z], with z decaying toward
% p.zfloor. This is identical to the model used for Experiments 1-5.

    p = pa_params();          % start from the full default set

    % turn OFF all post-scaffold features
    p.legacy_z  = 1;          % original single-z, fully-gated production
    p.s_dynamic = 0;          % algebraic salt signal (no tau_s lag)
    p.zfixed    = 0;          % no separate irreversible floor
    p.sodium    = 0;          % no dietary sodium offset
    p.zfloor    = 0.05;       % original capacity floor (z decays toward this)
end
