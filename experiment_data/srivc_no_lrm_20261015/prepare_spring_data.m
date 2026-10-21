function data=prepare_spring_data
% Extend the previous FRF data below 25 Hz, retaining its exact upper band.
base=fileparts(mfilename('fullpath')); out=fullfile(base,'spring_identification_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(fileparts(base),'auto_identify_motor1234_z0.mat'));
fs=1/median(diff(S.u_Motor1(:,1)));N=size(S.u_Motor1,1);U=zeros(N,4);Y=U;
for c=1:4
    U(:,c)=S.(sprintf('u_Motor%d',c))(:,2);
    Y(:,c)=S.(sprintf('pre_y_%d__0_',c-1))(:,2);
end
first=round(fs)+1;n=floor((N-first+1)/2);it=first:first+n-1;iv=first+n:first+2*n-1;
uscale=std(U(it,:),0,1);yscale=1e-6;fraw=(0:n-1)'*fs/n;
Ut=fft((U(it,:)-mean(U(it,:)))./uscale,n)/n;Yt=fft((Y(it,:)-mean(Y(it,:)))/yscale,n)/n;
Uv=fft((U(iv,:)-mean(U(iv,:)))./uscale,n)/n;Yv=fft((Y(iv,:)-mean(Y(iv,:)))/yscale,n)/n;
lowf=(.5:.1:24.9)'; windows=[15 30 45];low=struct('half_window',{},'train',{},'validation',{},'train_checks',{},'validation_checks',{});
for k=1:numel(windows)
    [ht,ct]=localFRF(fraw,Ut,Yt,lowf,uscale,yscale,windows(k));
    [hv,cv]=localFRF(fraw,Uv,Yv,lowf,uscale,yscale,windows(k));
    low(k)=struct('half_window',windows(k),'train',ht,'validation',hv,'train_checks',ct,'validation_checks',cv);
    sel=lowf>=3;
    fprintf('Low LRM half %d: repeatability 3-25 Hz %.4f%%, cond median/max %.3g/%.3g\n',windows(k), ...
        100*norm(ht(sel,:)-hv(sel,:),'fro')/norm(hv(sel,:),'fro'),median(ct.condition(:)),max(ct.condition(:)));
end
T=load(fullfile(base,'frf_preprocessing_comparison.mat'),'comparison');upper=T.comparison;
T=load(fullfile(fileparts(base),'analysis_20261015','auto_identify_motor1234_z0_spectral.mat'),'train','validation');
data=struct('low_f',lowf,'low',low,'upper',upper,'welchTrain',T.train,'welchValidation',T.validation, ...
    'fs',fs,'input_scale',uscale,'output_scale',yscale,'training_indices',[it(1) it(end)],'validation_indices',[iv(1) iv(end)], ...
    'source_file',fullfile(fileparts(base),'auto_identify_motor1234_z0.mat'), ...
    'note','Total closed-loop motor commands; d_z/d_Rx/d_Ry are identically zero, so not usable instruments. Validation reused for development.');
save(fullfile(out,'spring_frf_data.mat'),'data','-v7.3');
end

function [H,checks]=localFRF(f,U,Y,centers,uscale,yscale,half)
nu=4;ny=4;p=17;df=median(diff(f));H=complex(zeros(numel(centers),16));condition=zeros(numel(centers),4);ranks=condition;
for k=1:numel(centers)
    [~,center]=min(abs(f-centers(k)));first=max(2,min(center-half,numel(f)-2*half));ix=first:first+2*half;
    x=(f(ix)-centers(k))/(half*df);
    for j=1:ny
        X=complex(zeros(numel(ix),p));output=Y(ix,j);
        for power=0:2
            X(:,power*nu+(1:nu))=U(ix,:).*x.^power;X(:,nu*3+power+1)=x.^power;
        end
        X(:,end-1:end)=[-x.*output,-x.^2.*output];sc=max(vecnorm(X),realmin);Xs=X./sc;sv=svd(Xs,'econ');
        condition(k,j)=sv(1)/sv(end);ranks(k,j)=sum(sv>1e-12*sv(1));theta=lsqminnorm(Xs,output,1e-12)./sc.';
        for u=1:nu,H(k,j+(u-1)*ny)=theta(u)*yscale/uscale(u);end
    end
end
checks=struct('condition',condition,'rank',ranks,'parameter_count',p,'window_points',2*half+1,'window_span_hz',2*half*df);
assert(all(isfinite(H(:))) && all(ranks(:)==p),'Local FRF rank/finite check failed.');
end
