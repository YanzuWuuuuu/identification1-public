function report=verify_manuscript_derivatives()
% Independent checks of the equations added to the manuscript; no refitting.
base=fileparts(mfilename('fullpath'));addpath(base);
S=load(fullfile(base,'outputs_structured','distinct_N32','selected.mat'));
D=load(fullfile(base,'outputs','data_distinct.mat'));
x=S.selected.x;c=S.selected.chart;
[p,dm,dr]=clv_decode(x,c);
f=[22 100 286.43 344.58 488.27 705.67 798];ss=2i*pi*f;
Kf=freqresp(D.K,2*pi*f);exc=repmat(.1*eye(4),1,1,numel(f));
[~,J]=clv_mean(x,c,f,D.cfg.Ts,Kf,exc);
modalError=0;chartError=0;continuousError=[];sampledError=[];targetError=[];
for it=1:6
    v=sin((1:c.n)'*(it+.37));v=v/norm(v);h=1e-5;
    pp=clv_decode(x+h*v,c);pm=clv_decode(x-h*v,c);
    dM=sum(dm.*reshape(v(1:15),1,1,15),3);dR=zeros(4,4,c.K);
    for i=1:c.K
        ix=c.ro+14*(i-1)+(1:14);
        alpha=ones(4,1);alpha(c.order(2:4,i))=x(ix(1:3))+1i*x(ix(4:6));
        beta=x(ix(7:10))+1i*x(ix(11:14));
        Si=eye(4);Si=Si(:,c.order(2:4,i));
        JR=c.rs(i)*[kron(beta,Si),1i*kron(beta,Si),kron(eye(4),alpha),1i*kron(eye(4),alpha)];
        dR(:,:,i)=reshape(JR*v(ix),4,4);
        codeDeriv=sum(reshape(dr(:,:,i,:),16,14).*v(ix).',2);
        chartError=max(chartError,norm(JR*v(ix)-codeDeriv)/max(norm(codeDeriv),eps));
    end
    dl=p.lambda.*v(c.fo+(1:c.K))+p.omega.*(-p.zeta-1i*p.zeta.^2./sqrt(1-p.zeta.^2)).*v(c.zo+(1:c.K));
    de=reshape(c.ts*c.a_map*v(c.to+(1:7)),4,4);
    da=reshape(c.as*c.a_map*v(c.ao+(1:7)),4,4);db=p.b.*reshape(v(c.bo+(1:16)),4,4);
    dg=zeros(4,4,numel(f));
    for k=1:numel(f)
        z=ss(k);P=p.M/z^2;dP=dM/z^2;
        for i=1:c.K
            R=p.R(:,:,i);l=p.lambda(i);w=p.omega(i);zz=p.zeta(i);
            pair=R/(z-l)+conj(R)/(z-conj(l));
            rational=(2*real(R)*z+2*zz*w*real(R)-2*w*sqrt(1-zz^2)*imag(R))/(z^2+2*zz*w*z+w*w);
            modalError=max(modalError,norm(pair-rational,'fro')/max(norm(pair,'fro'),eps));
            P=P+pair;dP=dP+dR(:,:,i)/(z-l)+R*dl(i)/(z-l)^2+conj(dR(:,:,i))/(z-conj(l))+conj(R)*conj(dl(i))/(z-conj(l))^2;
        end
        den=1+p.a*z+p.b*z*z;ell=p.ty+p.tu.';
        dg(:,:,k)=exp(-ell*z)./den.*(dP-P.*(z*de+(z*da+z*z*db)./den));
    end
    fd=(continuous(pp,ss)-continuous(pm,ss))/(2*h);
    continuousError(it)=norm(fd(:)-dg(:))/norm(dg(:)); %#ok<AGROW>
    yp=clv_mean(x+h*v,c,f,D.cfg.Ts,Kf,exc);ym=clv_mean(x-h*v,c,f,D.cfg.Ts,Kf,exc);
    fd=(yp(:)-ym(:))/(2*h);sampledError(it)=norm(fd-J*v)/norm(J*v); %#ok<AGROW>
    [~,~,UU]=svd(p.M);nu=UU(:,4);[YY,~,~]=svd(p.M);ny=YY(:,4);
    dnu=-pinv(p.M)*dM*nu;dny=-pinv(p.M).'*dM.'*ny;
    dpiu=dnu*nu.'+nu*dnu.';dpiy=dny*ny.'+ny*dny.';
    analytic=[dpiy(:);dpiu(:)];
    for k=1:numel(ss)
        z=ss(k);den=1+p.a*z+p.b*z*z;deltaDen=z*da+z*z*db;ell=p.ty+p.tu.';
        ry=exp(-(ell(2:4,1)-ell(1,1))*z).*den(1,1)./den(2:4,1);
        ru=exp(-(ell(1,2:4)-ell(1,1))*z).*den(1,1)./den(1,2:4);
        dry=ry.*(-z*(de(2:4,1)-de(1,1))+deltaDen(1,1)/den(1,1)-deltaDen(2:4,1)./den(2:4,1));
        dru=ru.*(-z*(de(1,2:4)-de(1,1))+deltaDen(1,1)/den(1,1)-deltaDen(1,2:4)./den(1,2:4));
        dd=[dry;dru.'];analytic=[analytic;real(dd);imag(dd)]; %#ok<AGROW>
    end
    fd=(target(pp,ss)-target(pm,ss))/(2*h);targetError(it)=norm(fd-analytic)/norm(analytic); %#ok<AGROW>
end
% Test the derivative across an integer sample delay, not just within a cell.
ts=D.cfg.Ts;z=exp(ss*ts);hdelay=1e-9;
[H,dH]=clv_kernel(p.omega(1),p.zeta(1),p.a(1),p.b(1),ts,z,ts,true,true);
Hp=clv_kernel(p.omega(1),p.zeta(1),p.a(1),p.b(1),ts+hdelay,z,ts,true,false);
Hm=clv_kernel(p.omega(1),p.zeta(1),p.a(1),p.b(1),ts-hdelay,z,ts,true,false);
fd=(Hp-Hm)/(2*hdelay);dd=dH(:,:,5);
boundaryError=norm(fd-dd,'fro')/norm(dd,'fro');
% Check the profile least-squares identity with rank-revealing projections.
Jr=[real(J);imag(J)];nm=c.to;Jm=Jr(:,1:nm);Jc=Jr(:,nm+1:end);
dc=cos((1:size(Jc,2))');obs=sin((1:size(Jr,1))');rhs=obs-Jc*dc;
dmLS=lsqminnorm(Jm,rhs);Um=orth(Jm);left=norm(Jm*dmLS-rhs)^2;right=norm(rhs-Um*(Um'*rhs))^2;
report=struct('parameter_count',c.n,'directions',6,'modal_real_map_error',modalError,...
    'residue_kronecker_map_error',chartError,'continuous_jacobian_errors',continuousError,...
    'sampled_closed_jacobian_errors',sampledError,'calibration_target_jacobian_errors',targetError,...
    'integer_delay_boundary_error',boundaryError,'profile_projection_identity_error',abs(left-right)/max(left,eps),...
    'scope','Equation checks at saved estimate, no new fitting, no confidence-interval or convergence validation.');
assert(modalError<1e-10 && chartError<1e-10);
assert(max(continuousError)<1e-5 && max(sampledError)<1e-4 && max(targetError)<1e-5);
assert(boundaryError<1e-4 && report.profile_projection_identity_error<1e-8);
root=fileparts(fileparts(fileparts(base)));
fid=fopen(fullfile(root,'docs','research','辨识数学补全_导数核对_20261014.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(report,PrettyPrint=true));fclose(fid);
disp(jsonencode(report,PrettyPrint=true));
end

function G=continuous(p,s)
G=zeros(4,4,numel(s));ell=p.ty+p.tu.';
for k=1:numel(s)
    z=s(k);P=p.M/z^2;
    for i=1:numel(p.lambda),P=P+p.R(:,:,i)/(z-p.lambda(i))+conj(p.R(:,:,i))/(z-conj(p.lambda(i)));end
    G(:,:,k)=exp(-ell*z)./(1+p.a*z+p.b*z*z).*P;
end
end

function g=target(p,s)
[U,~,V]=svd(p.M);Py=U(:,4)*U(:,4).';Pu=V(:,4)*V(:,4).';g=[Py(:);Pu(:)];
ell=p.ty+p.tu.';
for k=1:numel(s)
    z=s(k);den=1+p.a*z+p.b*z*z;
    ry=exp(-(ell(2:4,1)-ell(1,1))*z).*den(1,1)./den(2:4,1);
    ru=exp(-(ell(1,2:4)-ell(1,1))*z).*den(1,1)./den(1,2:4);
    rr=[ry;ru.'];g=[g;real(rr);imag(rr)]; %#ok<AGROW>
end
end
