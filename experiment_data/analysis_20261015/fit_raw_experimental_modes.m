function fit_raw_experimental_modes(positions)
% Direct local frequency-domain output-error fits to raw U,Y FFTs.
% No spectral-matrix inversion, no nonparametric FRF in the objective.
if nargin<1,positions=[-.5 0 .5 1];end
out=fileparts(mfilename('fullpath'));base=fileparts(out);
windows={[240,315],[315,385],[455,535],[650,735]};
bounds={[260,300],[320,348;348,375],[475,510],[660,692;695,720]};
old=load(fullfile(out,'local_modal_models.mat'),'records');
allRecords=cell(numel(positions),1);rows={};groupRows={};
for ip=1:numel(positions)
    S=load(fullfile(base,sprintf('auto_identify_motor1234_z%g.mat',positions(ip))));
    t=S.u_Motor1(:,1);fs=1/median(diff(t));U=zeros(numel(t),4);Y=U;
    for k=1:4
        U(:,k)=S.(sprintf('u_Motor%d',k))(:,2);
        Y(:,k)=S.(sprintf('pre_y_%d__0_',k-1))(:,2);
    end
    start=round(fs)+1;n=floor((numel(t)-start+1)/2);nf=floor(n/2)+1;
    f=(0:nf-1)'*fs/n;
    Ut=fft(U(start:start+n-1,:)-mean(U(start:start+n-1,:)),n)/n;
    Yt=fft(Y(start:start+n-1,:)-mean(Y(start:start+n-1,:)),n)/n;
    Uv=fft(U(start+n:start+2*n-1,:)-mean(U(start+n:start+2*n-1,:)),n)/n;
    Yv=fft(Y(start+n:start+2*n-1,:)-mean(Y(start+n:start+2*n-1,:)),n)/n;
    oldIndex=find(cellfun(@(q)q.position_mm==positions(ip),old.records),1);
    groups=cell(4,1);
    for ig=1:4
        id=find(f>=windows{ig}(1)&f<=windows{ig}(2));
        % All DFT lines in each local band are retained.
        initial=old.records{oldIndex}.groups{ig}.train.frequency_hz;
        a=rawLocal(f(id),Ut(id,:),Yt(id,:),fs,windows{ig},bounds{ig},initial);
        b=rawLocal(f(id),Uv(id,:),Yv(id,:),fs,windows{ig},bounds{ig},a.frequency_hz);
        [Xg,Xtr]=rawDesign(f(id),Uv(id,:),fs,a.frequency_hz,a.damping,a.center_hz,a.scale_hz);
        forced=Xg*a.transfer_coefficients;
        % Only record-specific boundary/transient terms are re-estimated.
        wt=(f(id)/mean(f(id))).^2;
        A=[real(Xtr.*wt);imag(Xtr.*wt)];B=[real((Yv(id,:)-forced).*wt);imag((Yv(id,:)-forced).*wt)];
        sc=vecnorm(A);c=(A./sc)\B;c=c./sc.';
        predicted=forced+Xtr*c;
        validationError=norm((predicted-Yv(id,:)).*wt,'fro')/norm(Yv(id,:).*wt,'fro');
        forcedError=norm((forced-Yv(id,:)).*wt,'fro')/norm(Yv(id,:).*wt,'fro');
        groups{ig}=struct('train',a,'second_half',b,'f',f(id), ...
            'measured_validation',Yv(id,:),'predicted_validation',predicted, ...
            'validation_error',validationError,'forced_only_validation_error',forcedError);
        groupRows(end+1,:)={positions(ip),ig,100*a.relative_error,100*validationError,100*forcedError, ...
            a.design_condition,a.design_rank,a.design_columns,a.exitflag}; %#ok<AGROW>
        for j=1:numel(a.frequency_hz)
            rows(end+1,:)={positions(ip),ig,j,a.frequency_hz(j),100*a.damping(j), ...
                b.frequency_hz(j),100*b.damping(j),abs(a.frequency_hz(j)-b.frequency_hz(j)), ...
                100*validationError,a.residue_rank1_fraction(j),a.exitflag}; %#ok<AGROW>
        end
        fprintf('RAW z%g group%d f=%s damping%%=%s train=%.3f%% val=%.3f%% designCond=%.3g rank=%d/%d\n', ...
            positions(ip),ig,mat2str(a.frequency_hz.',6),mat2str(100*a.damping.',4), ...
            100*a.relative_error,100*validationError,a.design_condition,a.design_rank,a.design_columns);
    end
    allRecords{ip}=struct('position_mm',positions(ip),'groups',{groups},'sample_rate_hz',fs,'segment_seconds',n/fs);
end
modalTable=cell2table(rows,'VariableNames',{'z_mm','group','component','frequency_hz','damping_percent', ...
    'second_half_frequency_hz','second_half_damping_percent','split_frequency_difference_hz', ...
    'group_validation_error_percent','rank1_residue_fraction','exitflag'});
groupTable=cell2table(groupRows,'VariableNames',{'z_mm','group','training_error_percent','validation_error_percent', ...
    'forced_only_validation_error_percent','design_condition','design_rank','design_columns','exitflag'});
if numel(positions)==1,tag=sprintf('_z%g',positions);else,tag='';end
save(fullfile(out,['raw_local_models' tag '.mat']),'allRecords','modalTable','groupTable');
writetable(modalTable,fullfile(out,['raw_local_modes' tag '.csv']));
writetable(groupTable,fullfile(out,['raw_local_validation' tag '.csv']));
end

function fit=rawLocal(f,U,Y,fs,win,limits,start)
nm=size(limits,1);center=mean(f);scale=std(f);wt=(f/center).^2;
normalizer=norm(Y.*wt,'fro');
lb=log([limits(:,1);repmat(.00005,nm,1)]);ub=log([limits(:,2);repmat(.12,nm,1)]);
opts=optimoptions('lsqnonlin','Display','off','MaxIterations',100,'MaxFunctionEvaluations',1400, ...
    'FunctionTolerance',1e-10,'StepTolerance',1e-9,'OptimalityTolerance',1e-8);
best=Inf;
for damp=[.003,.025]
    x0=log([start(:);repmat(damp,nm,1)]);
    [xx,cost,~,flag,output]=lsqnonlin(@objective,min(max(x0,lb),ub),lb,ub,opts);
    if cost<best,best=cost;x=xx;exitflag=flag;iterations=output.iterations;end
end
[~,coef,pred,A]=objective(x);fn=exp(x(1:nm));damping=exp(x(nm+1:end));
sv=svd(A./vecnorm(A),'econ');nc=(2*nm+4)*4;rank1=zeros(nm,1);
for j=1:nm
    R=complex(zeros(4));
    for k=1:4
        idx=(k-1)*(2*nm+4)+2*j-1;
        R(:,k)=coef(idx,:).'+1i*coef(idx+1,:).';
    end
    sig=svd(R);rank1(j)=sig(1)^2/sum(sig.^2);
end
fit=struct('frequency_hz',fn,'damping',damping,'transfer_coefficients',coef(1:nc,:), ...
    'boundary_coefficients',coef(nc+1:end,:),'center_hz',center,'scale_hz',scale,'sample_rate_hz',fs, ...
    'window_hz',win,'relative_error',sqrt(best),'exitflag',exitflag,'iterations',iterations, ...
    'prediction',pred,'design_condition',sv(1)/sv(end),'design_rank',sum(sv>sv(1)*1e-10), ...
    'design_columns',size(A,2),'residue_rank1_fraction',rank1);
    function [r,coef,pred,A]=objective(x)
        [Xg,Xtr]=rawDesign(f,U,fs,exp(x(1:nm)),exp(x(nm+1:end)),center,scale);
        X=[Xg,Xtr];A=[real(X.*wt);imag(X.*wt)];B=[real(Y.*wt);imag(Y.*wt)];
        sc=vecnorm(A);coef=(A./sc)\B;coef=coef./sc.';
        pred=X*coef;r=(A*coef-B)/normalizer;r=r(:);
    end
end

function [Xg,Xtr]=rawDesign(f,U,fs,fn,damp,center,scale)
z=exp(2i*pi*f/fs);basis=complex(zeros(numel(f),2*numel(fn)+4));
for j=1:numel(fn)
    lambda=2*pi*fn(j)*(-damp(j)+1i*sqrt(1-damp(j)^2));p=exp(lambda/fs);
    plus=1./(z-p);minus=1./(z-conj(p));
    basis(:,2*j-1)=plus+minus;basis(:,2*j)=1i*(plus-minus);
end
v=(f-center)/scale;basis(:,end-3:end)=[ones(size(v)),1i*ones(size(v)),v,1i*v];
Xg=complex(zeros(numel(f),4*size(basis,2)));
for k=1:4,Xg(:,(k-1)*size(basis,2)+(1:size(basis,2)))=basis.*U(:,k);end
Xtr=basis;
end
