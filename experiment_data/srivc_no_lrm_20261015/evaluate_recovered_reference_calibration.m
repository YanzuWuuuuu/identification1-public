function e=evaluate_recovered_reference_calibration
% Check new spatial estimates on BOTH frontends; no check-half model choice.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'full_coordinate_calibration_20261016');
R=load(fullfile(out,'reference_identification.mat'),'ref');r=R.ref;
T=load(fullfile(out,'result.mat'),'result');prior=T.result;
T=load(fullfile(base,'structured_calibration_20261016','calibration_z0.mat'),'calib');c=T.calib;
f=r.f(:);nf=numel(f);nm=numel(r.reference_hz);target=prior.target_index;
P=prior.methods(1).Ty_static(1:3,:);Q=prior.methods(1).Tu_static(:,1:3);
ny0=prior.methods(1).Ty_static(4,:).';nu0=prior.methods(1).Tu_static(:,4);
oldny=prior.methods(3).Ty_static(4,:).';oldnu=prior.methods(3).Tu_static(:,4);
s0=spaces(r.direct_coefficients,r.a,r.reference_hz,target,ny0,nu0);
s1=spaces(r.moment_coefficients,r.a,r.reference_hz,target,ny0,nu0);
ny=s1.ny;nu=s1.nu;cf=s1.target.c;bf=s1.target.b;
ay=ny.'*cf;au=bf.'*nu;assert(abs(ay)>1e-8 && abs(au)>1e-8);
dP=-(P*cf)/ay*ny.';dQ=-nu*(bf.'*Q)/au;
ly=fraction(P,dP,ny,1.25,'y');lu=fraction(Q,dQ,nu,1.25,'u');
specs={'B0',P,Q,ny0,nu0,false; 'B4',P,Q,oldny,oldnu,true; ...
    'R0',P,Q,s0.ny,s0.nu,true; 'R1',P,Q,ny,nu,true; ...
    'R2',P+dP,Q+dQ,ny,nu,true; 'R3',P+ly*dP,Q+lu*dQ,ny,nu,true};
py=exp(-2i*pi*f*prior.extra_output_delay_seconds.');pu=exp(-2i*pi*f*prior.extra_input_delay_seconds);
data={reshape(c.H_validation.',4,4,[]),reshape(r.local{2}.FRF.',4,4,[])};
frontNames={'old_direct_FRF','recovered_reference_FRF'};methods=struct([]);gainRows={};
for m=1:size(specs,1)
    pm=specs{m,2};qm=specs{m,3};Y=[pm;specs{m,4}.'];U=[qm,specs{m,5}];
    methods(m).code=specs{m,1};methods(m).Ty=Y;methods(m).Tu=U; %#ok<AGROW>
    methods(m).Ty_physical=[pm.*prior.fixed_output_scales;specs{m,4}.'];
    methods(m).Tu_physical=[qm.*prior.fixed_input_scales,specs{m,5}];
    for source=1:2
        Z=complex(zeros(4,4,nf));
        for k=1:nf
            Cy=eye(4);Cu=eye(4);if specs{m,6},Cy=diag(py(k,:));Cu=diag(pu(k,:));end
            Z(:,:,k)=Y*Cy*data{source}(:,:,k)*Cu*U;
        end
        methods(m).responses{source}=Z;
    end
    gainRows(end+1,:)={specs{m,1},max(vecnorm(pm,2,2)./vecnorm(P,2,2)), ...
        max(vecnorm(qm,2,1)./vecnorm(Q,2,1)),norm(Y,2)/norm([P;ny.'],2), ...
        norm(U,2)/norm([Q,nu],2),cond(Y),cond(U)}; %#ok<AGROW>
    writematrix(methods(m).Ty_physical,fullfile(out,[specs{m,1} '_reference_Ty_output.csv']));
    writematrix(methods(m).Tu_physical,fullfile(out,[specs{m,1} '_reference_Tu_input.csv']));
end
checks=struct('estimated_support_output_basis_preserved',norm(dP*s1.C,'fro'), ...
    'estimated_support_input_basis_preserved',norm(s1.B*dQ,'fro'), ...
    'target_annihilation',[norm((P+dP)*cf),norm(bf.'*(Q+dQ))]);
checks.NN_invariant=0;
for source=1:2
    a=methods(4).responses{source}(4,4,:);
    for m=5:6,b=methods(m).responses{source}(4,4,:);checks.NN_invariant=max(checks.NN_invariant,norm(a(:)-b(:))/norm(a(:)));end
end
assert(max([checks.estimated_support_output_basis_preserved,checks.estimated_support_input_basis_preserved, ...
    checks.target_annihilation,checks.NN_invariant])<1e-9);
bands=[3 25;25 80;100 170;260 310;320 380;450 800;25 800];
labels={'3-25','25-80','100-170','260-310','320-380','450-800','25-800'};
blockNames={'logical_to_N','N_to_logical','N_to_N'};rows={};
for source=1:2
    for b=1:size(bands,1)
        ix=f>=bands(b,1) & f<=bands(b,2)+1e-9;fb=f(ix);
        for block=1:3
            v0=blockRms(fb,methods(1).responses{source},block,ix);
            v4=blockRms(fb,methods(2).responses{source},block,ix);
            for m=1:numel(methods)
                val=blockRms(fb,methods(m).responses{source},block,ix);
                rows(end+1,:)={frontNames{source},specs{m,1},labels{b},blockNames{block},val,val/v0,val/v4}; %#ok<AGROW>
            end
        end
    end
end
metrics=cell2table(rows,'VariableNames',{'frontend','method','band_hz','block','check_rms','over_B0','over_old_B4'});
gains=cell2table(gainRows,'VariableNames',{'method','max_logical_output_row_gain','max_logical_input_column_gain', ...
    'full_output_operator_vs_R1','full_input_operator_vs_R1','output_condition','input_condition'});

% Conditional check-half re-estimation with frozen poles, delays, weights.
% This diagnoses sensitivity and is NEVER used in any transform above.
q=1i*f./r.reference_hz.';den=1+q.*r.a(:,1).'+q.^2.*r.a(:,2).';
Phi=complex(zeros(nf,2*nm+1));for m=1:nm,Phi(:,2*m-1)=1./den(:,m);Phi(:,2*m)=q(:,m)./den(:,m);end
Phi(:,end)=1;np=size(Phi,2);coV=zeros(np,16);
for y=1:4
    X=complex(zeros(nf*4,np*4));rhs=reshape(squeeze(r.local{2}.Tyd(:,y,:)),[],1);
    for u=1:4
        puu=Phi.*exp(-2i*pi*f*r.delay(y,u));
        for d=1:4,X((d-1)*nf+(1:nf),(u-1)*np+(1:np))=puu.*r.local{2}.Tud(:,u,d);end
    end
    w=sqrt(repmat(r.moment_weights(:,y),4,1));Ar=X.*w;br=rhs.*w;
    Ar=[real(Ar);imag(Ar)];br=[real(br);imag(br)];sc=vecnorm(Ar);As=Ar./sc;
    assert(all(sc>0) && rank(As)==size(As,2));coef=(As\br)./sc.';
    for u=1:4,coV(:,y+(u-1)*4)=coef((u-1)*np+(1:np))*r.output_scale/r.input_scale(u);end
end
sv=spaces(coV,r.a,r.reference_hz,target,ny0,nu0);
sensitivity=struct('train_sigma4_over3',[s1.sy(4)/s1.sy(3),s1.su(4)/s1.su(3)], ...
    'check_sigma4_over3',[sv.sy(4)/sv.sy(3),sv.su(4)/sv.su(3)], ...
    'null_change_degrees',[angle(s1.ny,sv.ny),angle(s1.nu,sv.nu)], ...
    'target_change_degrees',[angle(s1.target.c,sv.target.c),angle(s1.target.b,sv.target.b)], ...
    'new_vs_old_null_degrees',[angle(ny,oldny),angle(nu,oldnu)]);
e=struct('f',f,'methods',methods,'metrics',metrics,'gains',gains,'source_direct',s0,'source_reference',s1, ...
    'source_check_reference',sv,'sensitivity',sensitivity,'target_projections',[ay au],'fractions',[ly lu], ...
    'checks',checks,'source_code',fileread([mfilename('fullpath') '.m']), ...
    'scope','Exploratory leading-three residue-space projection; no physical support truth or deployed feedback verification.');
save(fullfile(out,'reference_calibration_evaluation.mat'),'e','-v7.3');
writetable(metrics,fullfile(out,'reference_calibration_metrics.csv'));writetable(gains,fullfile(out,'reference_calibration_gains.csv'));
disp(sensitivity);disp(gains);disp(metrics(strcmp(metrics.frontend,'old_direct_FRF') & ...
    ismember(metrics.band_hz,{'25-80','260-310'}) & ~strcmp(metrics.block,'N_to_N'),:));
end
function s=spaces(co,a,refs,target,ny0,nu0)
Y=[];U=[];for m=1:3,R=residue(co,a,refs,m);Y=[Y,real(R),imag(R)];U=[U,real(R).',imag(R).'];end %#ok<AGROW>
[uy,sy]=svd(Y);[uu,su]=svd(U);ny=uy(:,end);nu=uu(:,end);
if ny.'*ny0<0,ny=-ny;end;if nu.'*nu0<0,nu=-nu;end
s=struct('C',uy(:,1:3),'B',uu(:,1:3).','ny',ny,'nu',nu,'sy',diag(sy),'su',diag(su),'target',realShape(residue(co,a,refs,target)));
end
function R=residue(co,a,refs,m)
p=roots([a(m,2),a(m,1),1]);p=p(imag(p)>0);
R=2*pi*refs(m)*(reshape(co(2*m-1,:),4,4)+p*reshape(co(2*m,:),4,4))/(a(m,1)+2*a(m,2)*p);
end
function t=realShape(R)
A=R/norm(R,'fro');g=linspace(0,pi,365);v=arrayfun(@(p)score(A,p),g);candidates=[0 pi];values=[v(1) v(end)];
ii=find(v(2:end-1)>=v(1:end-2)&v(2:end-1)>=v(3:end))+1;
for k=ii,[p,j]=fminbnd(@(x)-score(A,x),g(k-1),g(k+1),optimset('TolX',1e-11));candidates(end+1)=p;values(end+1)=-j;end %#ok<AGROW>
if v(1)>=v(2)&&v(1)>=v(end-1),[p,j]=fminbnd(@(x)-score(A,x),-g(2),g(2),optimset('TolX',1e-11));candidates(end+1)=p;values(end+1)=-j;end
[~,k]=max(values);[u,~,v]=svd(real(exp(-1i*candidates(k))*A));cv=u(:,1);bv=v(:,1);alpha=cv.'*R*bv;
t=struct('c',cv,'b',bv,'error',norm(R-alpha*cv*bv.','fro')/norm(R,'fro'));
end
function v=score(A,p)
s=svd(real(exp(-1i*p)*A));v=s(1)^2;
end
function l=fraction(A,D,n,cap,side)
if valid(A+D,A,n,cap,side),l=1;return;end
lo=0;hi=1;for k=1:70,m=(lo+hi)/2;if valid(A+m*D,A,n,cap,side),lo=m;else,hi=m;end;end;l=lo;
end
function ok=valid(A,B,n,cap,side)
if side=='y',g=vecnorm(A,2,2)./vecnorm(B,2,2);aa=[A;n.'];bb=[B;n.'];
else,g=vecnorm(A,2,1)./vecnorm(B,2,1);aa=[A,n];bb=[B,n];end
ok=max(g)<=cap && norm(A,2)<=cap*norm(B,2) && norm(aa,2)<=cap*norm(bb,2);
end
function a=angle(x,y)
a=acosd(min(1,abs(x.'*y)/norm(x)/norm(y)));
end
function val=blockRms(f,H,b,ix)
if b==1,z=H(4,1:3,ix);elseif b==2,z=H(1:3,4,ix);else,z=H(4,4,ix);end
x=reshape(z,[],nnz(ix)).';val=sqrt(trapz(f,sum(abs(x).^2,2))/(f(end)-f(1)));
end
