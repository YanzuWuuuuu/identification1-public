function campaign=prepare_position_spring_data
% Same split, local rational FRF settings and frequency grids at every z.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'position_spring_identification_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'spring_identification_20261016','spring_selected.mat'),'result');parent=S.result;
f=parent.f;low=f<25;positions=[-1 -.5 0 .5 1];datasets=cell(1,5);
for k=1:5
    tag=strrep(strrep(sprintf('%g',positions(k)),'-','m'),'.','p');
    path=fullfile(fileparts(base),sprintf('auto_identify_motor1234_z%g.mat',positions(k)));
    raw=load(path);fs=1/median(diff(raw.u_Motor1(:,1)));N=size(raw.u_Motor1,1);U=zeros(N,4);Y=U;
    for c=1:4
        U(:,c)=raw.(sprintf('u_Motor%d',c))(:,2);Y(:,c)=raw.(sprintf('pre_y_%d__0_',c-1))(:,2);
    end
    first=round(fs)+1;n=floor((N-first+1)/2);it=first:first+n-1;iv=first+n:first+2*n-1;
    uscale=std(U(it,:),0,1);yscale=1e-6;fr=(0:n-1)'*fs/n;
    Ut=fft((U(it,:)-mean(U(it,:)))./uscale,n)/n;Yt=fft((Y(it,:)-mean(Y(it,:)))/yscale,n)/n;
    Uv=fft((U(iv,:)-mean(U(iv,:)))./uscale,n)/n;Yv=fft((Y(iv,:)-mean(Y(iv,:)))/yscale,n)/n;
    % Upper local window uses the same truncated raw grid as the legacy z0
    % experiment. Low band is newly evaluated using the full positive DFT.
    ir=fr>=25&fr<=800;
    [Lt,CLt]=localFRF(fr,Ut,Yt,f(low),uscale,yscale,2);
    [Lv,CLv]=localFRF(fr,Uv,Yv,f(low),uscale,yscale,2);
    [Ht,CHt]=localFRF(fr(ir),Ut(ir,:),Yt(ir,:),f(~low),uscale,yscale,1);
    [Hv,CHv]=localFRF(fr(ir),Uv(ir,:),Yv(ir,:),f(~low),uscale,yscale,1);
    training=[Lt;Ht];validation=[Lv;Hv];df=[.1*ones(nnz(low),1);.5*ones(nnz(~low),1)];
    fl=max(1e-3*max(abs(Ht),[],1),eps);beta=ones(size(f));beta(low)=.03;
    W=(df.*beta)./max(abs(training),fl);W=W/median(W(:));
    time=raw.u_Motor1(:,1);tc=time-mean(time);slopes=(tc'*(U-mean(U,1)))/(tc'*tc);
    diag=struct('input_mean',mean(U,1),'input_std',std(U,0,1),'input_min',min(U,[],1), ...
        'input_max',max(U,[],1),'input_linear_slope_per_second',slopes,'sample_rate_hz',fs,'sample_count',N);
    warningText='';if positions(k)==-1,warningText='Motor 3 has an extreme drifting recorded command; diagnostic fit only, not a normal-position plant estimate.';end
    datasets{k}=struct('position_mm',positions(k),'tag',tag,'source_file',path,'f',f,'trainingFRF',training, ...
        'validationFRF',validation,'weights',W,'frequency_quadrature',df,'low_band_multiplier',.03, ...
        'input_scale',uscale,'output_scale',yscale,'raw_diagnostics',diag,'record_warning',warningText, ...
        'train_checks',combine(CLt,CHt));
    datasets{k}.train_checks=combine(CLt,CHt);datasets{k}.validation_checks=combine(CLv,CHv);
    datasets{k}.training_indices=[it(1) it(end)];datasets{k}.validation_indices=[iv(1) iv(end)];
    data=datasets{k};save(fullfile(out,['data_z' tag '.mat']),'data','-v7.3');
    fprintf('Prepared z=%g mm: low split %.3f%% upper split %.3f%% LRM condition median/max %.3g/%.3g\n',positions(k), ...
        100*norm(training(low,:)-validation(low,:),'fro')/norm(validation(low,:),'fro'), ...
        100*norm(training(~low,:)-validation(~low,:),'fro')/norm(validation(~low,:),'fro'), ...
        median(data.train_checks.condition(:)),max(data.train_checks.condition(:)));
end
campaign=struct('positions_mm',positions,'datasets',{datasets},'parent_file',fullfile(base,'spring_identification_20261016','spring_selected.mat'));
save(fullfile(out,'position_data.mat'),'campaign','-v7.3');
end

function [H,C]=localFRF(f,U,Y,centers,uscale,yscale,minimum)
half=45;nu=4;p=17;df=median(diff(f));H=complex(zeros(numel(centers),16));condition=zeros(numel(centers),4);ranks=condition;
for k=1:numel(centers)
    [~,center]=min(abs(f-centers(k)));first=max(minimum,min(center-half,numel(f)-2*half));ix=first:first+2*half;
    x=(f(ix)-centers(k))/(half*df);
    for j=1:4
        X=complex(zeros(numel(ix),p));output=Y(ix,j);
        for power=0:2
            X(:,power*nu+(1:nu))=U(ix,:).*x.^power;X(:,nu*3+power+1)=x.^power;
        end
        X(:,end-1:end)=[-x.*output,-x.^2.*output];sc=max(vecnorm(X),realmin);Xs=X./sc;sv=svd(Xs,'econ');
        condition(k,j)=sv(1)/sv(end);ranks(k,j)=sum(sv>1e-12*sv(1));theta=lsqminnorm(Xs,output,1e-12)./sc.';
        for u=1:4,H(k,j+(u-1)*4)=theta(u)*yscale/uscale(u);end
    end
end
C=struct('condition',condition,'rank',ranks,'parameter_count',p,'window_points',91,'window_span_hz',2*half*df);
assert(all(isfinite(H(:))) && all(ranks(:)==p),'Local FRF finite/rank check failed.');
end

function C=combine(L,H)
C=struct('condition',[L.condition;H.condition],'rank',[L.rank;H.rank],'parameter_count',17,'window_points',91,'window_span_hz',L.window_span_hz);
end
