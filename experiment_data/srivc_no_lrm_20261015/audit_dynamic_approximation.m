function bounds=audit_dynamic_approximation
% Unregularized real least-squares lower bound in each fixed basis family.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'dynamic_nullspace_20261016');
files=dir(fullfile(out,'*_p*.mat'));rows={};
for fi=1:numel(files)
    S=load(fullfile(files(fi).folder,files(fi).name),'result');r=S.result;
    f=r.f(:);ix=f>=25-1e-6 & f<=80+1e-6;f=f(ix);
    H=reshape(r.H_training(ix,:).',4,4,[]);n=numel(f);
    P=r.nominal.Ty./vecnorm(r.nominal.Ty,2,2);Q=r.nominal.Tu./vecnorm(r.nominal.Tu,2,1);
    Ey=[P;r.nominal.ny.'];Eu=[Q,r.nominal.nu];Z=complex(zeros(size(H)));
    for k=1:n,Z(:,:,k)=Ey*H(:,:,k)*Eu;end
    cfg=r.config;if ~isfield(cfg,'alignment_delay_seconds'),cfg.alignment_delay_seconds=0;end
    if ~isfield(cfg,'laguerre_order'),cfg.laguerre_order=0;end
    if ~isfield(cfg,'rich_basis'),cfg.rich_basis=false;end
    if cfg.pole_hz==0,phi=ones(n,1);else
        q=1i*f/cfg.pole_hz;phi=[1./(1+q).^2,q./(1+q).^2];
        if cfg.rich_basis,extra=phi(:,1)./(1+1i*f/20);phi=[phi,extra,extra./(1+1i*f/80)];end
        for ell=2:cfg.laguerre_order-1,phi(:,end+1)=(q-1).^ell./(q+1).^(ell+1);end %#ok<AGROW>
    end
    df=diff(f);w=[df(1);df(1:end-1)+df(2:end);df(end)]/2;nb=size(phi,2);
    for side=1:2
        M=complex(zeros(n*3,nb*3));
        for l=1:nb,for a=1:3
            if side==1,v=reshape(Z(a,1:3,:),3,[]).';else,v=reshape(Z(1:3,a,:),3,[]).';end
            v=v.*phi(:,l);M(:,a+3*(l-1))=v(:);
        end,end
        if side==1,v=reshape(Z(4,1:3,:),3,[]).';x=r.xy;block='logical_to_N';
        else,v=reshape(Z(1:3,4,:),3,[]).';x=r.xu;block='N_to_logical';end
        scale=sqrt(repmat(w,3,1)/sum(w.*sum(abs(v).^2,2)));
        v=v.*exp(-2i*pi*f*cfg.alignment_delay_seconds);v=v(:).*scale;M=M.*scale;
        A=[real(M);imag(M)];b=[real(v);imag(v)];xu=-(A\b);
        lower=norm(A*xu+b);actual=norm(A*x+b);
        check=abs(actual^2+cfg.ridge*sum(x.^2)-r.training_objective(side));
        assert(check<1e-9,'Independent objective disagrees with optimizer.');
        rows(end+1,:)={r.position_mm,r.tag,block,lower,actual,actual^2-lower^2,cond(A),check}; %#ok<AGROW>
    end
end
bounds=cell2table(rows,'VariableNames',{'position_mm','tag','block','fixed_basis_RMS_lower_bound','achieved_training_RMS','squared_cost_gap','real_design_condition','objective_check_error'});
save(fullfile(out,'approximation_bounds.mat'),'bounds');writetable(bounds,fullfile(out,'approximation_bounds.csv'));
disp(bounds(bounds.position_mm==0,:));
end
