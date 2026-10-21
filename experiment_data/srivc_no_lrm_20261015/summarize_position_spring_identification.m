function campaign=summarize_position_spring_identification
% Assemble all attempted positions, including the failed anomalous record.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'position_spring_identification_20261016');
positions=[-1 -.5 0 .5 1];results=cell(1,5);summaryRows={};channelRows={};modeRows={};delayRows={};
referenceBounds=[];referenceF=[];
for k=1:5
    tag=strrep(strrep(sprintf('%g',positions(k)),'-','m'),'.','p');
    S=load(fullfile(out,['model_z' tag '.mat']),'result');r=S.result;fit=r.selected.fit;f=r.f;
    if k==1,referenceBounds=r.common_bounds_hz;referenceF=f;end
    assert(isequal(referenceBounds,r.common_bounds_hz)&&isequal(referenceF,f),'Position settings are not comparable.');
    assert(isequal(size(r.trainingFRF),size(r.validationFRF),size(r.selected.H)) && all(isfinite(r.selected.H(:))));
    reconstructed=transfer(r,f);reconstruction=norm(reconstructed-r.selected.H,'fro')/norm(r.selected.H,'fro');
    ell=r.delay_matrix_seconds;additivity=max(abs(ell-ell(:,1)-ell(1,:)+ell(1,1)),[],'all');
    assert(reconstruction<1e-10 && additivity<1e-12);
    lo=f<25;hi=~lo;V=r.validationFRF;H=r.selected.H;
    vl=100*norm(H(lo,:)-V(lo,:),'fro')/norm(V(lo,:),'fro');vh=100*norm(H(hi,:)-V(hi,:),'fro')/norm(V(hi,:),'fro');
    repLow=100*norm(r.trainingFRF(lo,:)-V(lo,:),'fro')/norm(V(lo,:),'fro');
    repHigh=100*norm(r.trainingFRF(hi,:)-V(hi,:),'fro')/norm(V(hi,:),'fro');
    if isempty(r.selected.gradient),outerGradient=NaN;else,outerGradient=norm(r.selected.gradient,inf);end
    r.inner_converged=fit.converged;
    r.delay_optimization_complete=strcmp(r.outer_status,'completed') && r.optimizer_exitflag>0 && r.selected_matches_optimizer;
    if r.inner_converged&&r.delay_optimization_complete,r.status='converged';else,r.status='nonconverged';end
    r.verification=struct('reconstruction_relative_error',reconstruction,'delay_additivity_seconds',additivity, ...
        'inner_score',fit.score,'outer_gradient_inf',outerGradient);
    if ~isempty(r.record_warning),physicalStatus='anomalous_record';else,physicalStatus='nominal_record';end
    summaryRows(end+1,:)={positions(k),r.status,physicalStatus,fit.score,r.optimizer_exitflag,outerGradient, ...
        vl,vh,repLow,repHigh,median(r.train_checks.condition(:)),max(r.train_checks.condition(:)), ...
        min(r.train_checks.rank(:)),r.raw_diagnostics.input_linear_slope_per_second(3),reconstruction}; %#ok<AGROW>
    for y=1:4
        for u=1:4
            ch=y+4*(u-1);channelRows(end+1,:)={positions(k),y,u, ...
                100*norm(H(lo,ch)-V(lo,ch))/norm(V(lo,ch)),100*norm(H(hi,ch)-V(hi,ch))/norm(V(hi,ch)),r.status,physicalStatus}; %#ok<AGROW>
            delayRows(end+1,:)={positions(k),y,u,ell(y,u)*1e6,r.status,physicalStatus}; %#ok<AGROW>
        end
    end
    for m=1:numel(fit.frequency_hz)
        boundGap=min(abs(fit.frequency_hz(m)-r.common_bounds_hz(m,:)));
        modeRows(end+1,:)={positions(k),m,fit.frequency_hz(m),fit.damping(m),boundGap,r.status,physicalStatus}; %#ok<AGROW>
    end
    results{k}=r;result=r;save(fullfile(out,['model_z' tag '.mat']),'result','-v7.3');
end
summary=cell2table(summaryRows,'VariableNames',{'position_mm','numerical_status','record_status','inner_score','outer_exitflag','outer_gradient_inf', ...
    'validation_LRM_3_25_percent','validation_LRM_25_800_percent','FRF_split_3_25_percent','FRF_split_25_800_percent', ...
    'LRM_design_condition_median','LRM_design_condition_max','LRM_min_rank','motor3_linear_slope_per_s','model_reconstruction_relative_error'});
channels=cell2table(channelRows,'VariableNames',{'position_mm','output','input','validation_LRM_3_25_percent','validation_LRM_25_800_percent','numerical_status','record_status'});
modes=cell2table(modeRows,'VariableNames',{'position_mm','model_component_index','natural_frequency_hz','damping_ratio','nearest_frequency_bound_hz','numerical_status','record_status'});
delays=cell2table(delayRows,'VariableNames',{'position_mm','output','input','effective_total_delay_us','numerical_status','record_status'});
campaign=struct('results',{results},'positions_mm',positions,'summary',summary,'channel_metrics',channels,'modal_parameters',modes,'delay_table',delays, ...
    'common_frequency_grid_hz',referenceF,'common_bounds_hz',referenceBounds,'parent','spring_identification_20261016/spring_selected.mat', ...
    'question','Compare independently identified 3-800 Hz FRFs across the five recorded positions', ...
    'controls','Identical record split, LRM window, frequency grid, 15-component model family, inverse-magnitude weight rule, low-band coefficient 0.03, and seven-delay gauge/bounds', ...
    'limitations','z=-1 record abnormal and nonconverged; first/second halves reused exploratory validation; no physical mode matching or causal separation of position from run conditions', ...
    'execution','Local MATLAB R2025a, serial runs on i5-12600KF/32GB, no GPU, no raw file mutation', ...
    'next_action','Display the four converged normal-position fits together, retaining the anomalous attempt with explicit flags');
save(fullfile(out,'position_campaign.mat'),'campaign','-v7.3');
writetable(summary,fullfile(out,'position_identification_summary.csv'));
writetable(channels,fullfile(out,'position_channel_errors.csv'));
writetable(modes,fullfile(out,'position_modal_parameters.csv'));
writetable(delays,fullfile(out,'position_effective_delays.csv'));
disp(summary);
end

function H=transfer(r,f)
fit=r.selected.fit;nm=numel(fit.reference_hz);q=1i*f(:)./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';B=complex(zeros(numel(f),2*nm+1));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
B(:,end)=1;H=(B*fit.coefficients).*exp(-1i*2*pi*f(:).*r.delay_matrix_seconds(:).');
end
