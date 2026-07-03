function pa_mono_compare()
% PA_MONO_COMPARE  (clean version) Does the bistable salt bridge do independent
% work? Compares BISTABLE vs MONO-M1 (reversible) vs MONO-M2 (irreversible) on:
%   (1) DURABILITY  -- after a long pulse + withdrawal, settle to TRUE asymptote;
%                      remission = total capacity z stays LOW (disease cleared).
%   (2) RARITY      -- z_fixed threshold separating remitters/non-remitters.
%   (3) MIN DURATION-- shortest pulse that still durably remits.
%
% KEY FIX vs earlier: remission is judged on CAPACITY z=zfixed+z_p settling to
% a low stable asymptote (drug-off), NOT on a single aldosterone reading. This
% avoids being fooled by fast-output gating masking a capacity relapse.

    close all; here = fileparts(mfilename('fullpath'));
    opts = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',2.0);
    t_pre=50; t_treat=16*30; u_max=0.95;

    M = { 'BISTABLE',@pa_model,      pa_params();
          'MONO-M1', @pa_model_mono, pa_params_mono('reversible');
          'MONO-M2', @pa_model_mono, pa_params_mono('irreversible') };

    fprintf('\n%-10s | z_dis | z after pulse+withdrawal (TRUE asymptote) | durable capacity remission?\n','model');
    fprintf('%s\n',repmat('-',1,92));

    figure('Name','Clean monostable vs bistable','Color','w','Position',[40 40 1280 780]);
    for m=1:3
        name=M{m,1}; mdl=M{m,2}; pbase=M{m,3};

        % --- durability at remission-capable zfixed ---
        p=pbase; p.zfixed=0.05;
        xeq=settle(mdl,init_state(p,0.05,4.0,disease_zp(p)),p,opts);
        z_dis = p.zfixed + xeq(3);
        uft=@(t)u_max*(t>=t_pre && t<t_pre+t_treat);
        [t,X]=ode15s(@(t,x)mdl(t,x,p,uft),[0 t_pre+t_treat+200],xeq,opts);
        % from the post-treatment state, settle DRUG-OFF to the true asymptote
        xfin=settle(mdl,X(end,:).',p,opts);
        z_final = p.zfixed + xfin(3);
        durable = z_final < 0.5*z_dis;
        fprintf('%-10s | %5.2f | z_final=%.3f (z_p=%.3f) %s | %s\n', name, z_dis, ...
            z_final, xfin(3), repmat(' ',1,8), string(durable));

        % --- full trajectory for plotting (long follow) ---
        [tt,XX]=ode15s(@(t,x)mdl(t,x,p,uft),[0 t_pre+t_treat+4000],xeq,opts);
        subplot(3,2,2*m-1);
        plot(tt,XX(:,3)+p.zfixed,'LineWidth',1.8); hold on;     % total capacity z
        plot(tt,XX(:,2),'--','LineWidth',1.0);                   % aldosterone
        xline(t_pre,':'); xline(t_pre+t_treat,':');
        title(sprintf('%s: capacity z (solid) & a (dashed)',name));
        ylabel('z, a'); xlabel('days'); ylim([0 4.5]);
        yl=ylim; text(t_pre+t_treat+100, yl(2)*0.85, ...
            sprintf('z\\rightarrow%.2f',z_final),'FontSize',8);

        % --- rarity: z_fixed sweep (durable capacity remission) ---
        zfg=linspace(0,0.6,25); rem=false(size(zfg));
        for i=1:numel(zfg)
            q=pbase; q.zfixed=zfg(i);
            xe=settle(mdl,init_state(q,0.05,4.0,disease_zp(q)),q,opts);
            zd=q.zfixed+xe(3);
            uf=@(t)u_max*(t>=t_pre && t<t_pre+t_treat);
            [t2,X2]=ode15s(@(t,x)mdl(t,x,q,uf),[0 t_pre+t_treat+200],xe,opts);
            xf=settle(mdl,X2(end,:).',q,opts);
            rem(i) = (q.zfixed+xf(3)) < 0.5*zd;
        end
        thr=NaN; if any(rem)&&~all(rem), thr=zfg(find(rem,1,'last')); end
        subplot(3,2,2*m);
        bar(zfg,double(rem),1.0,'FaceColor',[0.4 0.7 0.4],'EdgeColor','none'); hold on;
        if ~isnan(thr), xline(thr,'r--','LineWidth',1.3); end
        title(sprintf('%s: durable remission vs z_{fixed}',name));
        xlabel('z_{fixed}'); ylabel('remits'); ylim([0 1.2]);

        % --- min duration ---
        dg=[15 30 60 120 240 480]; remd=false(size(dg));
        p=pbase; p.zfixed=0.05; xe=settle(mdl,init_state(p,0.05,4.0,disease_zp(p)),p,opts); zd=p.zfixed+xe(3);
        for j=1:numel(dg)
            uf=@(t)u_max*(t>=t_pre && t<t_pre+dg(j));
            [t3,X3]=ode15s(@(t,x)mdl(t,x,p,uf),[0 t_pre+dg(j)+200],xe,opts);
            xf=settle(mdl,X3(end,:).',p,opts);
            remd(j) = (p.zfixed+xf(3)) < 0.5*zd;
        end
        mind=NaN; if any(remd), mind=dg(find(remd,1,'first')); end
        fprintf('           rarity threshold z_fixed*=%s ; min duration=%s\n', ...
            tern(isnan(thr),'NONE',sprintf('%.3f',thr)), tern(isnan(mind),'never',sprintf('%dd',mind)));
    end

    sgtitle('Capacity-based remission: bistable durably clears disease; M1 regrows it; M2 clears by irreversibility');
    saveas(gcf,fullfile(here,'pa_mono_compare.png'));
    fprintf('\nSaved figure to pa_mono_compare.png\n');
end

function zp0=disease_zp(p), if isfield(p,'zp_set'), zp0=p.zp_set; else, zp0=max(p.z_disease-p.zfixed,0.1); end, end
function x0=init_state(p,r0,a0,zp0)
    if isfield(p,'s_dynamic')&&p.s_dynamic
        so=0; if isfield(p,'sodium'),so=p.sodium;end
        x0=[r0;a0;zp0;min(a0^p.ns/(p.Ks^p.ns+a0^p.ns)+so,0.999)];
    else, x0=[r0;a0;zp0]; end
end
function xe=settle(mdl,x0,p,opts), [~,X]=ode15s(@(t,x)mdl(t,x,p,@(tt)0),[0 12000],x0,opts); xe=X(end,:).'; end
function o=tern(c,a,b), if c,o=a;else,o=b;end, end
