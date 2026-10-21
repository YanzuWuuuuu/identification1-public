function evaluation=evaluate_dynamic_nullspace_revision
% Evaluate implemented stable filters; never fit/select on validation data.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'dynamic_nullspace_20261016');
S=load(fullfile(base,'band_calibration_20261016','band_calibration_campaign.mat'),'campaign');
campaign=S.campaign;bands=[3 25;25 80;80 170;100 170;260 310;320 380;450 800;25 800];
rows={};checks={};responses={};
for p=1:4
    c=campaign.records{p};f=c.f(:);H=reshape(c.H_validation.',4,4,[]);
    P=c.nominal.Ty./vecnorm(c.nominal.Ty,2,2);Q=c.nominal.Tu./vecnorm(c.nominal.Tu,2,1);
    Ey=[P;c.nominal.ny.'];Eu=[Q,c.nominal.nu];Z0=transform(H,Ey,Eu);
    for m=1:numel(campaign.results{p}.methods)
        old=campaign.results{p}.methods(m);
        if ~ismember(old.code,{'B0','B1','B4','D0','D2'}),continue;end
        addMetrics(old.H_validation,Z0,c.position_mm,old.code,'own_position',f);
        responses{end+1}=struct('position_mm',c.position_mm,'method',old.code,'scope','own_position','f',f,'Z',old.H_validation); %#ok<AGROW>
    end
    files=dir(fullfile(out,[c.tag '_p*.mat']));
    for m=1:numel(files)
        T=load(fullfile(files(m).folder,files(m).name),'result');r=T.result;
        code=erase(r.tag,[c.tag '_']);
        [Z,check]=applyFilters(H,r);
        addMetrics(Z,Z0,c.position_mm,code,'own_position',f,delayOf(r));
        checks(end+1,:)={c.position_mm,code,'own_position',r.accepted,r.filters.output_hinf,r.filters.input_hinf,check(1),check(2),check(3),check(4),check(5)}; %#ok<AGROW>
        responses{end+1}=struct('position_mm',c.position_mm,'method',code,'scope','own_position','f',f,'Z',Z); %#ok<AGROW>
    end
    if p~=2
        files=dir(fullfile(out,'z0_p*.mat'));
        for m=1:numel(files)
            T=load(fullfile(files(m).folder,files(m).name),'result');r=T.result;
            assert(isequal(f,r.f),'Cross-position frequency grids differ.');
            code=erase(r.tag,'z0_');[Z,check]=applyFilters(H,r);
            addMetrics(Z,Z0,c.position_mm,code,'z0_transfer',f,delayOf(r));
            checks(end+1,:)={c.position_mm,code,'z0_transfer',r.accepted,r.filters.output_hinf,r.filters.input_hinf,check(1),check(2),check(3),check(4),check(5)}; %#ok<AGROW>
            responses{end+1}=struct('position_mm',c.position_mm,'method',code,'scope','z0_transfer','f',f,'Z',Z); %#ok<AGROW>
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'position_mm','method','scope','fmin','fmax','block','rms','over_B0','complex_change_over_B0','transported_complex_change'});
checkTable=cell2table(checks,'VariableNames',{'position_mm','method','scope','numerically_accepted','output_hinf','input_hinf','realization_error','anchor_error','A_relative_error','output_min_norm','input_min_norm'});
evaluation=struct('metrics',metrics,'checks',checkTable,'responses',{responses},'bands',bands, ...
    'note','Second-half LRM diagnostics already inspected during earlier work; not blind repeats, not isolated rigid/flexible truth.');
save(fullfile(out,'dynamic_evaluation.mat'),'evaluation','-v7.3');
writetable(metrics,fullfile(out,'dynamic_metrics.csv'));writetable(checkTable,fullfile(out,'implementation_checks.csv'));
ix=metrics.position_mm==0 & metrics.fmin==25 & metrics.fmax==80 & ~strcmp(metrics.block,'N_to_N');
disp(metrics(ix,{'method','block','over_B0'}));
ix=metrics.position_mm==0 & metrics.fmin>=260 & strcmp(metrics.block,'N_to_N');
disp(metrics(ix,{'method','fmin','fmax','complex_change_over_B0'}));

    function addMetrics(Z,Z0,pos,code,scope,f,alignmentDelay)
        if nargin<7,alignmentDelay=0;end
        for bi=1:size(bands,1)
            band=bands(bi,:);ix=f>=band(1)-1e-6 & f<=band(2)+1e-6;w=quadrature(f(ix));
            for bb=1:3
                if bb==1,x=reshape(Z(4,1:3,ix),3,[]).';x0=reshape(Z0(4,1:3,ix),3,[]).';block='logical_to_N';
                elseif bb==2,x=reshape(Z(1:3,4,ix),3,[]).';x0=reshape(Z0(1:3,4,ix),3,[]).';block='N_to_logical';
                else,x=reshape(Z(4,4,ix),[],1);x0=reshape(Z0(4,4,ix),[],1);block='N_to_N';end
                rms=@(v)sqrt(sum(w.*sum(abs(v).^2,2))/sum(w));b0=rms(x0);
                if bb==3,transport=exp(-4i*pi*f(ix)*alignmentDelay);else,transport=exp(-2i*pi*f(ix)*alignmentDelay);end
                rows(end+1,:)={pos,code,scope,band(1),band(2),block,rms(x),rms(x)/b0,rms(x-x0)/b0,rms(x-x0.*transport)/b0}; %#ok<AGROW>
            end
        end
    end
end
function [Z,ch]=applyFilters(H,r)
f=r.f(:);P=r.nominal.Ty./vecnorm(r.nominal.Ty,2,2);Q=r.nominal.Tu./vecnorm(r.nominal.Tu,2,1);
ny=r.nominal.ny(:);nu=r.nominal.nu(:);Ey=[P;ny.'];Eu=[Q,nu];Z0=transform(H,Ey,Eu);
fy=freqresp(r.filters.fy,2*pi*f);fu=freqresp(r.filters.fu,2*pi*f);
alignment=exp(-2i*pi*f*delayOf(r));Ny=complex(zeros(1,4,numel(f)));Nu=complex(zeros(4,1,numel(f)));
for k=1:numel(f),Ny(:,:,k)=alignment(k)*ny.'+fy(:,:,k)*P;Nu(:,:,k)=alignment(k)*nu+Q*fu(:,:,k);end
if delayOf(r)==0
    assert(max(abs(Ny-freqresp(r.filters.Ny,2*pi*f)),[],'all')<1e-8);
    assert(max(abs(Nu-freqresp(r.filters.Nu,2*pi*f)),[],'all')<1e-8);
end
phi=independentBasis(f,r.config);py=phi*reshape(r.xy,3,[]).';pu=phi*reshape(r.xu,3,[]).';
realErr=max([max(abs(reshape(fy,3,[]).'-py),[],'all'),max(abs(reshape(fu,3,[]).'-pu),[],'all')]);
assert(realErr<1e-8,'Actual stable realization differs from fitted frequency basis.');
Z=complex(zeros(size(H)));anchorErr=0;yn=zeros(numel(f),1);un=yn;
for k=1:numel(f)
    L=[eye(3),zeros(3,1);fy(:,:,k),alignment(k)];R=[eye(3),fu(:,:,k);zeros(1,3),alignment(k)];
    Z(:,:,k)=L*Z0(:,:,k)*R;
    direct=[P;Ny(:,:,k)]*H(:,:,k)*[Q,Nu(:,:,k)];
    assert(norm(Z(:,:,k)-direct,'fro')<1e-10*max(norm(Z(:,:,k),'fro'),realmin),'Triangular/physical paths disagree.');
    anchorErr=max(anchorErr,max(abs([Ny(:,:,k)*ny,nu.'*Nu(:,:,k)]-alignment(k))));
    yn(k)=norm(Ny(:,:,k));un(k)=norm(Nu(:,:,k));
end
aerr=norm(reshape(Z(1:3,1:3,:)-Z0(1:3,1:3,:),[],1))/norm(reshape(Z0(1:3,1:3,:),[],1));
assert(anchorErr<1e-10 && aerr<1e-14 && min(yn)>1-1e-10 && min(un)>1-1e-10);
ch=[realErr,anchorErr,aerr,min(yn),min(un)];
end
function Z=transform(H,Ey,Eu)
Z=complex(zeros(size(H)));for k=1:size(H,3),Z(:,:,k)=Ey*H(:,:,k)*Eu;end
end
function phi=independentBasis(f,cfg)
if cfg.pole_hz==0,phi=ones(numel(f),1);return;end
if ~isfield(cfg,'rich_basis'),cfg.rich_basis=false;end
if ~isfield(cfg,'laguerre_order'),cfg.laguerre_order=0;end
q=1i*f/cfg.pole_hz;phi=[1./(1+q).^2,q./(1+q).^2];
if cfg.rich_basis
    p2=phi(:,1)./(1+1i*f/20);phi=[phi,p2,p2./(1+1i*f/80)];
end
for ell=2:cfg.laguerre_order-1,phi(:,end+1)=((q-1)./(q+1)).^ell./(1+q);end %#ok<AGROW>
end
function w=quadrature(f)
df=diff(f);w=[df(1);df(1:end-1)+df(2:end);df(end)]/2;
end
function d=delayOf(r)
d=0;if isfield(r.config,'alignment_delay_seconds'),d=r.config.alignment_delay_seconds;end
end
