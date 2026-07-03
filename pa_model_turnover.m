function dx = pa_model_turnover(t, x, p, u_of_t)
% PA_MODEL_TURNOVER  Two-compartment DIRECT-TURNOVER comparator (the "fair"
% no-feedback competitor requested in review).
%
% This is the strongest no-bistability alternative to the salt-bridge model.
% Autonomous capacity is split into TWO patient-specific compartments:
%
%   z_R : RESISTANT capacity  -- drug-insensitive, irreversible (mutation-
%         locked). Constant. Analogous to z_fixed in the main model.
%   z_S : SENSITIVE capacity  -- driven down DIRECTLY by the drug (CYP11B2
%         inhibition -> adrenal turnover/involution), one-way (no regrowth),
%         NO feedback, NO salt dependence.
%
% Total autonomous capacity z = z_R + z_S. A "patient" is defined by the total
% disease capacity z_tot and the RESISTANT FRACTION f_R = z_R/z_tot:
%       z_R      = f_R * z_tot           (set in params, constant)
%       z_S(0)   = (1 - f_R) * z_tot     (initial condition)
%
% Drug action on sensitive compartment (irreversible turnover):
%       dz_S = eps * ( -g_u * u * z_S )         (one-way; u=0 -> dz_S=0)
%
% After a sufficiently long full-dose course, z_S -> 0, so the post-withdrawal
% total capacity -> z_R. Whether the patient "remits" then depends ONLY on
% whether z_R is low enough that aldosterone normalizes. There is NO positive
% feedback and NO bistability: for any fixed u the system has a single
% equilibrium, and the post-withdrawal state is a CONTINUOUS function of z_R.
%
% The point of the comparison (pa_turnover_compare.m): does remission appear as
% a SHARP threshold (a discontinuity, as in the bistable model's saddle-node)
% or as a SMOOTH gradient in the patient-intrinsic resistant fraction? A smooth
% boundary here, versus a discontinuous one in the bistable model, is a genuine
% qualitative discriminator between the two mechanisms.
%
%   States [r; a; z_S; s]   (z_R is a constant parameter p.zR, not a state)
%
% Fast plant (r,a) and salt state s are IDENTICAL to the main model so the
% comparison is like-for-like. Aldosterone output is LINEAR in total capacity
% (no phi(s) gating) so that the comparison is scored on capacity without the
% fast-loop output gating confounding it.

    r   = x(1);
    a   = x(2);
    z_S = x(3);

    sodium = 0; if isfield(p,'sodium'), sodium = p.sodium; end
    S_target = a^p.ns/(p.Ks^p.ns + a^p.ns) + sodium;
    S_target = min(max(S_target,0),0.999);
    if isfield(p,'s_dynamic') && p.s_dynamic && numel(x) >= 4
        s = x(4);
    else
        s = S_target;
    end

    zR = p.zR;                          % resistant (constant) capacity
    z  = zR + max(z_S,0);               % total autonomous capacity
    u  = u_of_t(t);

    % Output linear in total capacity (no phi(s) gating): no fast-loop bistability.
    auto_prod  = p.alpha * z;
    renin_prod = p.beta * r^p.pr/(p.Kr^p.pr + r^p.pr);

    dr = p.r0*(1-s) - p.kr*r;
    da = (1-u)*(renin_prod + auto_prod) - p.ka*a;

    % Sensitive compartment: direct, irreversible drug-driven turnover. No
    % feedback, no salt dependence, no regrowth (one-way).
    dz_S = p.eps * ( -p.gu * u * z_S );

    if isfield(p,'s_dynamic') && p.s_dynamic && numel(x) >= 4
        ds = (S_target - s)/p.tau_s;
        dx = [dr; da; dz_S; ds];
    else
        dx = [dr; da; dz_S];
    end
end
