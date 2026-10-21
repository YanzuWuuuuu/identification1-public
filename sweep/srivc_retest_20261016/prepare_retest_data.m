function campaign=prepare_retest_data
% New records, original data untouched. Quadratic local rational FRFs.
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
source=fullfile(root,'sweep','sweep');
names={'sweep_z_minus0_6.mat','sweep_z_minus0_3.mat','sweep_z_0.mat','sweep_z_plus0_5.mat','sweep_z_plus0_10.mat'};
positions=[-.6,-.3,0,.5,1];tags={'m0p6','m0p3','0','0p5','1'};
f=[(3:.1:24.9)';(25:.5:800)'];df=[.1*ones(220,1);.5*ones(1551,1)];
assert(numel(df)==numel(f));records=cell(1,5);rows={};
for k=1:5
    path=fullfile(source,names{k});raw=load(path);t=raw.SaturationU_0__0_(:,1);N=numel(t);U=zeros(N,4);Y=U;
    for j=0:3
        a=raw.(sprintf('SaturationU_0__%d_',j));b=raw.(sprintf('y_sensor_0__%d_',j));
        assert(isequal(size(a),[N 2])&&isequal(size(b),[N 2]));
        assert(max(abs(a(:,1)-t))<1e-10 && max(abs(b(:,1)-t))<1e-10);
        U(:,j+1)=a(:,2);Y(:,j+1)=b(:,2);
    end
    rawU=U;inputRecordFactor=(10/6)/12.9;
    % User confirmed this extra factor was present in the new motor records.
    % Undo the multiplication before FRF estimation: H_corrected=k*H_recorded.
    U=U/inputRecordFactor;
    fs=round(1/median(diff(t)));timeQuantization=max(abs(t-t(1)-(0:N-1)'/fs));
    % Exported timestamps are rounded to 0.1 microsecond; no missing samples.
    assert(timeQuantization<2e-7 && all(isfinite([U(:);Y(:)])));
    first=round(fs)+1;n=floor((N-first+1)/2);it=first:first+n-1;iv=first+n:first+2*n-1;
    us=std(U(it,:),0,1);ys=1e-6;fr=(0:n-1)'*fs/n;assert(all(us>0));
    Ut=fft((U(it,:)-mean(U(it,:)))./us,n)/n;Yt=fft((Y(it,:)-mean(Y(it,:)))/ys,n)/n;
    Uv=fft((U(iv,:)-mean(U(iv,:)))./us,n)/n;Yv=fft((Y(iv,:)-mean(Y(iv,:)))/ys,n)/n;
    low=f<25;hi=fr>=25 & fr<=800;
    [ht1,ct1]=localfrf(fr,Ut,Yt,f(low),us,ys,2);
    [hv1,cv1]=localfrf(fr,Uv,Yv,f(low),us,ys,2);
    [ht2,ct2]=localfrf(fr(hi),Ut(hi,:),Yt(hi,:),f(~low),us,ys,1);
    [hv2,cv2]=localfrf(fr(hi),Uv(hi,:),Yv(hi,:),f(~low),us,ys,1);
    H=[ht1;ht2];V=[hv1;hv2];floorval=max(1e-3*max(abs(ht2),[],1),eps);beta=ones(size(f));beta(low)=.03;
    W=(df.*beta)./max(abs(H),floorval);W=W/median(W(:));
    tc=t-mean(t);slopes=(tc.'*(U-mean(U)))/(tc.'*tc);
    diagnostics=struct('fs',fs,'N',N,'duration',t(end)-t(1),'time_start',t(1),'time_quantization_s',timeQuantization, ...
        'input_record_multiplier_removed',inputRecordFactor,'raw_mean_U',mean(rawU),'raw_std_U',std(rawU),'mean_U',mean(U),'std_U',std(U), ...
        'min_U',min(U),'max_U',max(U),'mean_Y',mean(Y),'std_Y',std(Y),'input_slope',slopes, ...
        'std_train_U',std(U(it,:)),'std_check_U',std(U(iv,:)),'std_train_Y',std(Y(it,:)),'std_check_Y',std(Y(iv,:)));
    r=struct('position_mm',positions(k),'tag',tags{k},'f',f,'trainingFRF',H,'validationFRF',V, ...
        'weights',W,'frequency_quadrature',df,'input_scale',us,'output_scale',ys,'source_file',path, ...
        'training_indices',[it(1) it(end)],'validation_indices',[iv(1) iv(end)],'diagnostics',diagnostics, ...
        'train_condition',[ct1.condition;ct2.condition],'check_condition',[cv1.condition;cv2.condition], ...
        'min_rank',min([ct1.rank(:);ct2.rank(:);cv1.rank(:);cv2.rank(:)]),'window_span_hz',90*fs/n, ...
        'record_warning','','method','Closed-loop total U/Y, local rational FRF; no independent external reference in new files');
    records{k}=r;data=r;save(fullfile(out,['data_z' tags{k} '.mat']),'data','-v7.3');
    loerr=100*norm(H(low,:)-V(low,:),'fro')/norm(V(low,:),'fro');hierr=100*norm(H(~low,:)-V(~low,:),'fro')/norm(V(~low,:),'fro');
    rows(end+1,:)={positions(k),names{k},N,fs,t(end)-t(1),loerr,hierr,median(r.train_condition(:)),max(r.train_condition(:)),r.min_rank,r.window_span_hz}; %#ok<AGROW>
    fprintf('Prepared z=%g: split low %.3f%%, high %.3f%%, local cond %.3g/%.3g\n',positions(k),loerr,hierr,median(r.train_condition(:)),max(r.train_condition(:)));
end
audit=cell2table(rows,'VariableNames',{'position_mm','file','sample_count','fs','duration_s','split_3_25_pct','split_25_800_pct','condition_median','condition_max','minimum_rank','window_span_hz'});
campaign=struct('records',{records},'positions_mm',positions,'raw_audit',audit,'source_code',fileread([mfilename('fullpath') '.m']),'matlab_version',version);
save(fullfile(out,'data_campaign.mat'),'campaign','-v7.3');writetable(audit,fullfile(out,'raw_and_frf_audit.csv'));disp(audit);
end
function [H,C]=localfrf(f,U,Y,centers,us,ys,minimum)
half=45;p=17;df=median(diff(f));H=complex(zeros(numel(centers),16));condi=zeros(numel(centers),4);ranks=condi;
for k=1:numel(centers)
    [~,center]=min(abs(f-centers(k)));first=max(minimum,min(center-half,numel(f)-2*half));ix=first:first+2*half;x=(f(ix)-centers(k))/(half*df);
    for j=1:4
        X=complex(zeros(numel(ix),p));out=Y(ix,j);
        for pow=0:2,X(:,pow*4+(1:4))=U(ix,:).*x.^pow;X(:,13+pow)=x.^pow;end
        X(:,16:17)=[-x.*out,-x.^2.*out];sc=max(vecnorm(X),realmin);Xs=X./sc;sv=svd(Xs,'econ');
        condi(k,j)=sv(1)/sv(end);ranks(k,j)=sum(sv>1e-12*sv(1));theta=lsqminnorm(Xs,out,1e-12)./sc.';
        for u=1:4,H(k,j+(u-1)*4)=theta(u)*ys/us(u);end
    end
end
assert(all(isfinite(H(:)))&&all(ranks(:)==p),'Local rank/finiteness failed.');C=struct('condition',condi,'rank',ranks);
end
