function checks=verify_fir_export()
% Verify saved coefficients, both implementations, and independent FRF arrays.
folder=fileparts(mfilename('fullpath'));s=load(fullfile(folder,'selected_arrays.mat'));
sources={'nonparam','srivc'};rows={};rng(20261017);
for k=1:2
    name=sources{k};r=load_dynamic_fir_trial(name);f=s.f(:);
    Fu=squeeze(freqresp(r.Fu,2*pi*f)).';Fy=squeeze(freqresp(r.Fy,2*pi*f)).';
    targetU=s.(['Fu_' name]);targetY=s.(['Fy_' name]);
    errU=norm(Fu-targetU,'fro')/max(norm(targetU,'fro'),1e-14);
    errY=norm(Fy-targetY,'fro')/max(norm(targetY,'fro'),1e-14);
    q=squeeze(freqresp(r.q,2*pi*f));wanted=exp(-2i*pi*f*r.added_delay_samples_per_side*r.Ts);
    errQ=max(abs(q(:)-wanted));
    H=s.(['P_' name]);P=s.Pstatic;multError=0;
    for j=1:numel(f)
        U=eye(4);Y=eye(4);U(1:3,4)=-Fu(j,:).';U(4,4)=q(j);
        Y(4,1:3)=-Fy(j,:);Y(4,4)=q(j);
        one=U; %#ok<NASGU>
        got=Y*squeeze(P(j,:,:))*U;target=squeeze(H(j,:,:));
        multError=max(multError,norm(got-target,'fro')/max(norm(target,'fro'),realmin));
    end
    n=1500;t=(0:n-1)'*r.Ts;un=randn(n,1);yr=randn(n,3);
    gotU=lsim(r.Fu,un,t);wantU=zeros(n,3);
    for j=1:3,wantU(:,j)=filter(r.hu(:,j),1,un);end
    gotY=lsim(r.Fy,yr,t);wantY=zeros(n,1);
    for j=1:3,wantY=wantY+filter(r.hy(:,j),1,yr(:,j));end
    timeU=norm(gotU-wantU,'fro')/max(norm(wantU,'fro'),1e-14);
    timeY=norm(gotY-wantY)/max(norm(wantY),1e-14);
    assert(max([errU,errY,errQ,multError,timeU,timeY])<1e-9,'Export mismatch.');
    assert(isstable(r.Ud)&&isstable(r.Yd),'FIR maps unstable.');
    rows(end+1,:)={name,errU,errY,errQ,multError,timeU,timeY}; %#ok<AGROW>
end
checks=cell2table(rows,'VariableNames',{'source','Fu_frequency_relerror','Fy_frequency_relerror', ...
    'delay_abs_error','full_product_relerror','Fu_time_relerror','Fy_time_relerror'});
writetable(checks,fullfile(folder,'implementation_checks.csv'));disp(checks);
end
