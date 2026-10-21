function report=summarize_delay_srivc
base=fileparts(mfilename('fullpath'));out=fullfile(base,'delay_identification_20261016');
Q=load(fullfile(out,'common_delay_srivc.mat'),'result');C=Q.result;
Q=load(fullfile(out,'channel_delay_srivc.mat'),'result');A=Q.result;
Q=load(fullfile(base,'srivc_z0_mid_fixed_full.mat'),'result');raw=Q.result;
models=struct('tag',{'no_delay','common_delay','additive_channel_delay'}, ...
    'fit',{C.baseline,C.selected.fit,A.selected.fit},'H',{C.baselineH,C.selected.H,A.selected.H}, ...
    'delay',{zeros(4),C.delay_matrix_seconds,A.delay_matrix_seconds});
f=C.f;W=C.weights;rows={};localRows={};channelRows={};
for k=1:numel(models)
    m=models(k);h=m.H;iv=f>=200&f<=735;
    errors=100*sqrt(sum(abs(h-C.welchValidation).^2,1)./sum(abs(C.welchValidation).^2,1));
    for ch=1:16,channelRows(end+1,:)={m.tag,mod(ch-1,4)+1,ceil(ch/4),errors(ch)};end %#ok<AGROW>
    ib=f>=180&f<=260;
    for ch=[1 6]
        fb=f(ib);[~,loc]=min(abs(C.trainingFRF(ib,ch)));center=fb(loc);
        masks={ib,ib&abs(f-center)<=5,ib&abs(f-center)>5};names={'180_260','center_plusminus5','remaining_flanks'};
        for b=1:3
            mask=masks{b};ref=C.lrmValidation(mask,ch);pred=h(mask,ch);
            localRows(end+1,:)={m.tag,ch,names{b},100*norm(pred-ref)/norm(ref)}; %#ok<AGROW>
        end
        [minimum,index]=min(abs(h(ib,ch)));
        localRows(end+1,:)={m.tag,ch,'minimum_db',20*log10(minimum)}; %#ok<AGROW>
        localRows(end+1,:)={m.tag,ch,'minimum_hz',fb(index)}; %#ok<AGROW>
    end
    % Validate in the original record: transfer frozen, same boundary family.
    hf=evaluateTransfer(m.fit,raw.validation.f,m.delay);
    forced=complex(zeros(size(raw.validation.Y)));
    for u=1:4
        forced=forced+hf(:,(u-1)*4+(1:4)).*(raw.validation.U(:,u)*raw.input_scale(u)/raw.output_scale);
    end
    % Boundary terms correspond to the undelayed rational model per output.
    % For delayed systems this is only an auxiliary flexible nuisance fit;
    % report it explicitly, never call it a full zero-state prediction test.
    [B,~,~,~,~,bc]=raw_srivc_design(raw.validation.f,raw.validation.U,m.fit.a,m.fit.reference_hz,true);
    L=B(:,bc);w=(raw.validation.f/100).^2;X=[real(L.*w);imag(L.*w)];
    y=(raw.validation.Y-forced).*w;sc=max(vecnorm(X),realmin);coeff=((X./sc)\[real(y);imag(y)])./sc.';
    pred=forced+L*coeff;den=norm(raw.validation.Y.*w,'fro');
    rows(end+1,:)={m.tag,100*norm((h-C.trainingFRF).*sqrt(W),'fro')/norm(C.trainingFRF.*sqrt(W),'fro'), ...
        100*norm((h-C.lrmValidation).*sqrt(W),'fro')/norm(C.lrmValidation.*sqrt(W),'fro'), ...
        100*norm(h-C.welchValidation,'fro')/norm(C.welchValidation,'fro'), ...
        100*norm(h(iv,:)-C.welchValidation(iv,:),'fro')/norm(C.welchValidation(iv,:),'fro'), ...
        100*norm(h-C.lrmValidation,'fro')/norm(C.lrmValidation,'fro'), ...
        100*norm((pred-raw.validation.Y).*w,'fro')/den, ...
        100*norm((forced-raw.validation.Y).*w,'fro')/den, min(errors),median(errors),max(errors)}; %#ok<AGROW>
end
summary=cell2table(rows,'VariableNames',{'model','weighted_train_error_percent','weighted_validation_lrm_percent', ...
    'welch_frf_error_25_800_percent','welch_frf_error_200_735_percent','lrm_frf_error_25_800_percent', ...
    'auxiliary_boundary_refitted_raw_output_percent','forced_raw_output_percent', ...
    'minimum_channel_error_percent','median_channel_error_percent','maximum_channel_error_percent'});
local=cell2table(localRows,'VariableNames',{'model','channel','metric','value'});
channels=cell2table(channelRows,'VariableNames',{'model','output','input','welch_relative_complex_error_percent'});
modeTable=table((1:numel(C.baseline.frequency_hz)).',C.baseline.frequency_hz,C.selected.fit.frequency_hz,A.selected.fit.frequency_hz, ...
    100*C.baseline.damping,100*C.selected.fit.damping,100*A.selected.fit.damping, ...
    'VariableNames',{'component','no_delay_hz','common_delay_hz','channel_delay_hz','no_delay_damping_percent','common_delay_damping_percent','channel_delay_damping_percent'});
report=struct('summary',summary,'local',local,'channel_errors',channels,'modes',modeTable, ...
    'delay_matrix_microseconds',1e6*A.delay_matrix_seconds,'outer_exitflag',A.optimizer_exitflag, ...
    'inner_score',A.selected.fit.score,'outer_gradient',A.selected.gradient, ...
    'raw_output_metric','Auxiliary validation only: record-specific rational boundary terms refitted; no exact delayed sampled boundary derivation');
save(fullfile(out,'delay_validation_metrics.mat'),'report');
writetable(summary,fullfile(out,'delay_validation_summary.csv'));writetable(local,fullfile(out,'delay_validation_notches.csv'));
writetable(channels,fullfile(out,'delay_validation_channels.csv'));writetable(modeTable,fullfile(out,'delay_modal_parameters.csv'));
writematrix(1e6*A.delay_matrix_seconds,fullfile(out,'effective_delay_matrix_us.csv'));
disp(summary);disp(local);disp(modeTable);
end

function H=evaluateTransfer(p,f,delay)
q=1i*f(:)./p.reference_hz(:).';den=1+q.*p.a(:,1).'+q.^2.*p.a(:,2).';
B=complex(zeros(numel(f),2*size(p.a,1)+1));
for m=1:size(p.a,1),B(:,2*m-1)=1./den(:,m);B(:,2*m)=q(:,m)./den(:,m);end
B(:,end)=1;H=(B*p.coefficients).*exp(-1i*2*pi*f(:).*delay(:).');
end
