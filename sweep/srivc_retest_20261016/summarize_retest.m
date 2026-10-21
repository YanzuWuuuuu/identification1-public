function campaign=summarize_retest
% Audit and summarize independent SRIVC fits of the five new sweep records.
% Does not refit, select a model, update delays, or overwrite model_z*.mat.
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
oldbase=fullfile(root,'experiment_data','srivc_no_lrm_20261015');
S=load(fullfile(oldbase,'position_spring_identification_20261016','position_campaign.mat'),'campaign');
old=S.campaign.results;
S=load(fullfile(oldbase,'static_decoupling_20261016','decoupled_campaign.mat'),'pack');
pack=S.pack;assert(strcmp(pack.views(2).tag,'geometric'));
Ty=pack.views(2).Ty;Tu=pack.views(2).Tu;
assert(isequal(size(Ty),[4 4])&&isequal(size(Tu),[4 4])&&rank(Ty)==4&&rank(Tu)==4);
positions=[-.6 -.3 0 .5 1];tags={'m0p6','m0p3','0','0p5','1'};
bands=[3 25;25 80;100 170;200 320;320 400;450 800;25 800];
bandNames={'3-25','25-80','100-170','200-320','320-400','450-800','25-800'};
records=cell(1,5);nullRecords=cell(1,5);comparisons=cell(1,3);
summaryRows={};modeRows={};bandRows={};channelRows={};matchRows={};changeRows={};changeChannelRows={};delayRows={};scaleRows={};scaleChannelRows={};
for k=1:numel(tags)
    S=load(fullfile(out,['model_z' tags{k} '.mat']),'result');r=S.result;fit=r.selected.fit;f=r.f(:);
    assert(abs(r.position_mm-positions(k))<1e-9 && all(diff(f)>0));
    assert(isequal(size(r.trainingFRF),[numel(f) 16]) && isequal(size(r.validationFRF),[numel(f) 16]));
    assert(all(isfinite(r.trainingFRF(:)))&&all(isfinite(r.validationFRF(:))));
    assert(strcmp(fit.method,'frequency_domain_additive_SRIVC_FRF'),'Expected actual additive SRIVC.');
    H=transfer(r,f);reconstruction=norm(H-r.selected.H,'fro')/max(norm(H,'fro'),realmin);
    ell=r.delay_matrix_seconds;additivity=max(abs(ell-ell(:,1)-ell(1,:)+ell(1,1)),[],'all');
    assert(reconstruction<1e-9 && additivity<1e-12,'Saved model or additive delays failed.');
    assert(r.diagnostics.time_quantization_s<2e-7,'Unexpected timestamp deviation from uniform sample grid.');
    assert(isfield(r.diagnostics,'input_record_multiplier_removed') && ...
        abs(r.diagnostics.input_record_multiplier_removed-(10/6)/12.9)<1e-12, ...
        'This summary requires models reidentified after the user-confirmed motor-record factor correction.');
    nm=numel(fit.frequency_hz);assert(isequal(size(fit.bounds_hz),[nm 2]));
    q=r.diagnostics;delayStatus=fieldOr(r,'delay_status','unspecified');
    summaryRows(end+1,:)={r.position_mm,r.status,fit.method,nm,fit.converged,fit.score,fit.score_tolerance,fit.reason,fit.iterations, ...
        r.selected.start_index,fit.weighted_cost,delayStatus,q.fs,q.N,q.duration,1e6*q.time_quantization_s, ...
        median(r.train_condition(:)),max(r.train_condition(:)),median(r.check_condition(:)),max(r.check_condition(:)), ...
        r.min_rank,r.window_span_hz,reconstruction,additivity,q.input_record_multiplier_removed}; %#ok<AGROW>
    for m=1:nm
        fn=fit.frequency_hz(m);gap=min(fn-fit.bounds_hz(m,1),fit.bounds_hz(m,2)-fn);
        modeRows(end+1,:)={r.position_mm,m,fn,fit.damping(m),100*fit.damping(m),fit.bounds_hz(m,1),fit.bounds_hz(m,2), ...
            gap,fit.score,fit.score_tolerance,fit.converged,r.status,delayStatus}; %#ok<AGROW>
    end
    LH=transform(H,Ty,Tu);LT=transform(r.trainingFRF,Ty,Tu);LV=transform(r.validationFRF,Ty,Tu);
    groups={1:16,[1 2 3 5 6 7 9 10 11],16,[4 8 12],[13 14 15]};
    names={'physical_4x4','logical_3x3','null_to_null','logical_to_null','null_to_logical'};
    for b=1:size(bands,1)
        mask=bandMask(f,bands(b,:));fb=f(mask);
        for g=1:numel(groups)
            if g==1,hh=H(mask,:);tt=r.trainingFRF(mask,:);vv=r.validationFRF(mask,:);
            else,ii=groups{g};hh=LH(mask,ii);tt=LT(mask,ii);vv=LV(mask,ii);end
            rt=complexRMS(fb,tt);rv=complexRMS(fb,vv);et=complexRMS(fb,hh-tt);ev=complexRMS(fb,hh-vv);es=complexRMS(fb,tt-vv);
            bandRows(end+1,:)={r.position_mm,names{g},bandNames{b},numel(fb),100*et/max(rt,realmin), ...
                100*ev/max(rv,realmin),100*es/max(rv,realmin),rt,rv,et,ev,es,r.status,delayStatus}; %#ok<AGROW>
        end
        for view=1:2
            if view==1,hh=H;tt=r.trainingFRF;vv=r.validationFRF;label='physical';else,hh=LH;tt=LT;vv=LV;label='geometric_full';end
            for u=1:4
                for y=1:4
                    c=y+4*(u-1);rv=complexRMS(fb,vv(mask,c));ev=complexRMS(fb,hh(mask,c)-vv(mask,c));es=complexRMS(fb,tt(mask,c)-vv(mask,c));
                    channelRows(end+1,:)={r.position_mm,label,y,u,bandNames{b},100*ev/max(rv,realmin),100*es/max(rv,realmin),rv,ev,r.status}; %#ok<AGROW>
                end
            end
        end
    end
    for u=1:4
        for y=1:4,delayRows(end+1,:)={r.position_mm,y,u,1e6*ell(y,u),delayStatus};end %#ok<AGROW>
    end
    records{k}=r;
    nullRecords{k}=struct('position_mm',r.position_mm,'f',f,'model_full_geometric',LH,'training_full_geometric',LT,'check_full_geometric',LV, ...
        'model_NN',LH(:,16),'training_NN',LT(:,16),'check_NN',LV(:,16), ...
        'model_logical_to_N',LH(:,[4 8 12]),'check_logical_to_N',LV(:,[4 8 12]), ...
        'model_N_to_logical',LH(:,[13 14 15]),'check_N_to_logical',LV(:,[13 14 15]), ...
        'status',r.status,'delay_status',delayStatus);
end

% Only exact shared positions are paired. Negative positions have changed.
pairPositions=[0 .5 1];
for k=1:numel(pairPositions)
    z=pairPositions(k);on=find(cellfun(@(r)abs(r.position_mm-z)<1e-9,old));nn=find(abs(positions-z)<1e-9);
    assert(isscalar(on)&&isscalar(nn));o=old{on};n=records{nn};
    matchRows=[matchRows;matchModes(o,n)]; %#ok<AGROW>
    % Compare only genuinely common measured bins. In particular, the new
    % 25 Hz bin has no old upper-band partner; do not interpolate across the
    % independently fitted low/high LRM boundary at 24.9--25.5 Hz.
    % The old FFT sampling-rate calculation leaves about 3.3e-11 relative
    % roundoff in its bin labels. One microhertz is still 1e5 times smaller
    % than the finest 0.1 Hz grid step and cannot merge neighboring bins.
    f=n.f(:);oldFrequency=o.f(:);[mask,oldIndex]=ismembertol(f,oldFrequency,1e-6,'DataScale',1);f=f(mask);oldIndex=oldIndex(mask);
    assert(all(oldIndex>0)&&max(abs(f-oldFrequency(oldIndex)))<1e-6&&numel(unique(oldIndex))==numel(oldIndex));
    O=o.validationFRF(oldIndex,:);N=n.validationFRF(mask,:);
    OT=o.trainingFRF(oldIndex,:);NT=n.trainingFRF(mask,:);
    OH=transfer(o,f);NH=transfer(n,f);
    LO=transform(O,Ty,Tu);LN=transform(N,Ty,Tu);LOT=transform(OT,Ty,Tu);LNT=transform(NT,Ty,Tu);
    LOH=transform(OH,Ty,Tu);LNH=transform(NH,Ty,Tu);
    assert(all(isfinite(O(:)))&&all(isfinite(LO(:))));
    % The user-confirmed extra command factor has already been removed upstream.
    % Check remaining common scaling in 40--80 Hz as a separate diagnostic only.
    ss=f>=40&f<80;oa=O(ss,:);na=N(ss,:);
    commonScale=sum(conj(oa(:)).*na(:))/sum(abs(oa(:)).^2);
    assert(isfinite(commonScale)&&abs(commonScale)>realmin);
    scaleRows(end+1,:)={z,40,80,nnz(ss),real(commonScale),imag(commonScale),abs(commonScale),angle(commonScale)*180/pi, ...
        20*log10(abs(commonScale)),100*norm(na/commonScale-oa,'fro')/max(norm(oa,'fro'),realmin),o.status,n.status}; %#ok<AGROW>
    for u=1:4
        for y=1:4
            ch=y+4*(u-1);ratio=na(:,ch)./oa(:,ch);assert(all(isfinite(ratio)));
            phase=unwrap(angle(ratio))*180/pi;phase=mod(median(phase)+180,360)-180;
            scaleChannelRows(end+1,:)={z,y,u,40,80,median(abs(ratio)),20*log10(median(abs(ratio))),phase, ...
                min(abs(ratio)),max(abs(ratio)),o.status,n.status}; %#ok<AGROW>
        end
    end
    for b=1:size(bands,1)
        ii=bandMask(f,bands(b,:));fb=f(ii);
        for view=1:2
            if view==1,a=O;bnew=N;at=OT;bt=NT;ah=OH;bh=NH;label='physical';else,a=LO;bnew=LN;at=LOT;bt=LNT;ah=LOH;bh=LNH;label='geometric_full';end
            ar=complexRMS(fb,a(ii,:));br=complexRMS(fb,bnew(ii,:));dr=complexRMS(fb,bnew(ii,:)-a(ii,:));
            splita=complexRMS(fb,at(ii,:)-a(ii,:));splitb=complexRMS(fb,bt(ii,:)-bnew(ii,:));
            modelChange=complexRMS(fb,bh(ii,:)-ah(ii,:));modelOld=complexRMS(fb,ah(ii,:));
            changeRows(end+1,:)={z,label,bandNames{b},100*dr/max(ar,realmin),ar,br,dr,20*log10(max(br,realmin)/max(ar,realmin)), ...
                100*splita/max(ar,realmin),100*splitb/max(br,realmin),100*modelChange/max(modelOld,realmin), ...
                100*complexRMS(fb,bnew(ii,:)/commonScale-a(ii,:))/max(ar,realmin), ...
                100*complexRMS(fb,ah(ii,:)-a(ii,:))/max(ar,realmin), ...
                100*complexRMS(fb,bh(ii,:)-bnew(ii,:))/max(br,realmin),o.status,n.status}; %#ok<AGROW>
            for u=1:4
                for y=1:4
                    c=y+4*(u-1);arc=complexRMS(fb,a(ii,c));brc=complexRMS(fb,bnew(ii,c));drc=complexRMS(fb,bnew(ii,c)-a(ii,c));
                    changeChannelRows(end+1,:)={z,label,y,u,bandNames{b},100*drc/max(arc,realmin), ...
                        20*log10(max(brc,realmin)/max(arc,realmin)),arc,brc,drc, ...
                        100*complexRMS(fb,bnew(ii,c)/commonScale-a(ii,c))/max(arc,realmin),o.status,n.status}; %#ok<AGROW>
                end
            end
        end
    end
    comparisons{k}=struct('position_mm',z,'f',f,'old_check_physical',O,'new_check_physical',N,'old_model_physical',OH,'new_model_physical',NH, ...
        'old_check_geometric',LO,'new_check_geometric',LN,'old_model_geometric',LOH,'new_model_geometric',LNH, ...
        'common_scale_40_80_diagnostic',commonScale,'scale_diagnostic','Remaining-scale diagnostic after user-confirmed U_corrected=U_recorded/((10/6)/12.9); no further normalization of primary FRFs, models, or figures', ...
        'old_status',o.status,'new_status',n.status);
end

summary=cell2table(summaryRows,'VariableNames',{'position_mm','numerical_status','method','component_count','inner_converged','inner_score','score_tolerance','stop_reason', ...
    'iterations','selected_start','weighted_training_cost','delay_status','fs_hz','sample_count','duration_s','time_grid_rounding_us', ...
    'training_condition_median','training_condition_max','check_condition_median','check_condition_max','minimum_LRM_rank','LRM_window_span_hz','reconstruction_error','delay_additivity_error_s','input_record_multiplier_removed'});
modes=cell2table(modeRows,'VariableNames',{'position_mm','component_index','frequency_hz','damping_ratio','damping_percent','frequency_lower_bound_hz','frequency_upper_bound_hz', ...
    'nearest_frequency_bound_hz','inner_score','score_tolerance','inner_converged','numerical_status','delay_status'});
metrics=cell2table(bandRows,'VariableNames',{'position_mm','response_group','band_hz','frequency_count','training_complex_RMS_pct','check_complex_RMS_pct','split_complex_RMS_pct', ...
    'training_FRF_RMS','check_FRF_RMS','training_error_RMS','check_error_RMS','split_error_RMS','numerical_status','delay_status'});
channels=cell2table(channelRows,'VariableNames',{'position_mm','coordinates','output_index','input_index','band_hz','check_complex_RMS_pct','split_complex_RMS_pct','check_FRF_RMS','check_error_RMS','numerical_status'});
matching=cell2table(matchRows,'VariableNames',{'position_mm','match_status','old_component','new_component','old_frequency_hz','new_frequency_hz','frequency_change_pct', ...
    'old_damping_ratio','new_damping_ratio','damping_change_percentage_points','damping_relative_change_pct','residue_MAC','delay_removed_residue_MAC', ...
    'frequency_relative_distance','greedy_match_cost','old_gate_candidate_count','new_gate_candidate_count','old_status','new_status'});
changes=cell2table(changeRows,'VariableNames',{'position_mm','coordinates','band_hz','new_old_complex_FRF_change_pct','old_check_FRF_RMS','new_check_FRF_RMS','difference_RMS', ...
    'new_old_RMS_ratio_dB','old_split_complex_RMS_pct','new_split_complex_RMS_pct','model_new_old_complex_RMS_change_pct','common_scale_removed_shape_error_pct', ...
    'old_model_check_complex_RMS_pct','new_model_check_complex_RMS_pct','old_status','new_status'});
changeChannels=cell2table(changeChannelRows,'VariableNames',{'position_mm','coordinates','output_index','input_index','band_hz','new_old_complex_FRF_change_pct', ...
    'new_old_RMS_ratio_dB','old_check_FRF_RMS','new_check_FRF_RMS','difference_RMS','common_scale_removed_shape_error_pct','old_status','new_status'});
delays=cell2table(delayRows,'VariableNames',{'position_mm','output_index','input_index','effective_delay_us','delay_status'});
scales=cell2table(scaleRows,'VariableNames',{'position_mm','frequency_low_hz','frequency_high_excluded_hz','frequency_count','common_scale_real','common_scale_imag', ...
    'common_scale_magnitude','common_scale_phase_deg','common_scale_dB','scale_removed_40_80_shape_error_pct','old_status','new_status'});
scaleChannels=cell2table(scaleChannelRows,'VariableNames',{'position_mm','output_index','input_index','frequency_low_hz','frequency_high_excluded_hz', ...
    'median_magnitude_ratio','median_magnitude_ratio_dB','median_unwrapped_ratio_phase_deg','minimum_magnitude_ratio','maximum_magnitude_ratio','old_status','new_status'});

campaign=struct('records',{records},'positions_mm',positions,'summary',summary,'modal_parameters',modes,'band_metrics',metrics,'channel_metrics',channels, ...
    'matched_old_new_modes',matching,'old_new_FRF_changes',changes,'old_new_channel_changes',changeChannels,'delay_table',delays, ...
    'common_scale_diagnostic',scales,'channel_scale_diagnostic',scaleChannels, ...
    'static_output_map',Ty,'static_input_map',Tu,'axis_labels',{{'Rx','Ry','Z','N'}},'bands_hz',bands, ...
    'units','Recorded sensor units / user-corrected actuator commands: U_corrected=U_recorded/((10/6)/12.9). The known extra motor-record multiplier has been removed before identification. Static logical scaling uses prior nominal geometry; no independent absolute mechanical-unit calibration is claimed.', ...
    'RMS_definition','sqrt(trapz(f,sum(abs(H).^2,2))/(f(end)-f(1))); upper band edge excluded except 800 Hz. Five-new-record metrics retain all new bins; old/new pair metrics use only measured-frequency intersections at 1e-6 Hz absolute tolerance (old sampling-rate roundoff is about 3.3e-11 relative), with no interpolation across the LRM split.', ...
    'mode_matching','Greedy one-to-one ascending df/f_old/0.1+(1-MAC)/0.3, gated df/f_old<0.1 and complex observed-pole-residue MAC>0.7; unmatched modes retained; no physical labels guaranteed', ...
    'scale_diagnostic','After the user-confirmed motor factor correction, one remaining common complex LS scale is diagnosed on physical check FRFs at 40<=f<80 Hz. The divided-new shape residual is a separate diagnostic only; there is no further normalization of the primary FRFs, modal parameters, or figures. The earlier approximately 7.74-fold difference was an explained recording-factor mismatch, not evidence of a screw-induced mechanical gain change.', ...
    'validation_scope','Reused two halves of each closed-loop record; total U/Y FRFs not certified unbiased; same-position pre/post checks do not isolate torque from all other changes', ...
    'timestamp_handling','Original exported timestamps rounded to 0.1 microsecond; prepare_retest_data uses integer sampling rate and verified uniform sample count; no timestamp resampling', ...
    'source_code',fileread([mfilename('fullpath') '.m']),'matlab_version',version);
save(fullfile(out,'campaign.mat'),'campaign','-v7.3');
writetable(summary,fullfile(out,'identification_summary.csv'));writetable(modes,fullfile(out,'modal_parameters.csv'));
writetable(metrics,fullfile(out,'band_metrics.csv'));writetable(channels,fullfile(out,'channel_metrics.csv'));
writetable(matching,fullfile(out,'matched_old_new_modes.csv'));writetable(changes,fullfile(out,'old_new_FRF_changes.csv'));
writetable(changeChannels,fullfile(out,'old_new_channel_changes.csv'));writetable(delays,fullfile(out,'effective_delays.csv'));
writetable(scales,fullfile(out,'old_new_common_scale_diagnostic.csv'));writetable(scaleChannels,fullfile(out,'old_new_channel_scale_diagnostic.csv'));
writematrix(Ty,fullfile(out,'geometric_Ty_full.csv'));writematrix(Tu,fullfile(out,'geometric_Tu_full.csv'));
writematrix(Ty(1:3,:),fullfile(out,'geometric_Ty_logical.csv'));writematrix(Tu(:,1:3),fullfile(out,'geometric_Tu_logical.csv'));
writematrix(Ty(4,:).',fullfile(out,'output_null_vector.csv'));writematrix(Tu(:,4),fullfile(out,'input_null_vector.csv'));
nullspaces=struct('records',{nullRecords},'same_position_comparisons',{comparisons},'Ty',Ty,'Tu',Tu,'axis_labels',{{'Rx','Ry','Z','N'}}, ...
    'matrix_layout','Nfrequency x 16, column major: output varies fastest; input logical->N columns 4,8,12; N->logical columns 13,14,15', ...
    'operation','Fixed static geometry only, no delay compensation; model delays retained','units',campaign.units);
save(fullfile(out,'nullspace_responses.mat'),'nullspaces','-v7.3');disp(summary);
end

function rows=matchModes(o,n)
fo=o.selected.fit;fn=n.selected.fit;[Ro,R0o]=residues(o);[Rn,R0n]=residues(n);
no=numel(fo.frequency_hz);nn=numel(fn.frequency_hz);d=abs(fo.frequency_hz(:)-fn.frequency_hz(:).')./fo.frequency_hz(:);
mac=zeros(no,nn);mac0=mac;
for i=1:no
    for j=1:nn,mac(i,j)=MAC(Ro(:,i),Rn(:,j));mac0(i,j)=MAC(R0o(:,i),R0n(:,j));end
end
gate=d<.1 & mac>.7;cost=d/.1+(1-mac)/.3;cost(~gate)=Inf;
[sorted,indices]=sort(cost(:));usedOld=false(no,1);usedNew=false(nn,1);rows={};
for p=1:numel(indices)
    if ~isfinite(sorted(p)),break;end
    [i,j]=ind2sub(size(cost),indices(p));if usedOld(i)||usedNew(j),continue;end
    usedOld(i)=true;usedNew(j)=true;
    rows(end+1,:)=matchRow(o,n,i,j,'matched',mac(i,j),mac0(i,j),d(i,j),cost(i,j),sum(gate(i,:)),sum(gate(:,j))); %#ok<AGROW>
end
for i=find(~usedOld).',rows(end+1,:)=matchRow(o,n,i,NaN,'unmatched_old',NaN,NaN,NaN,NaN,sum(gate(i,:)),0);end %#ok<AGROW>
for j=find(~usedNew).',rows(end+1,:)=matchRow(o,n,NaN,j,'unmatched_new',NaN,NaN,NaN,NaN,0,sum(gate(:,j)));end %#ok<AGROW>
end

function row=matchRow(o,n,i,j,status,mac,mac0,d,cost,co,cn)
of=NaN;oz=NaN;nf=NaN;nz=NaN;
if isfinite(i),of=o.selected.fit.frequency_hz(i);oz=o.selected.fit.damping(i);end
if isfinite(j),nf=n.selected.fit.frequency_hz(j);nz=n.selected.fit.damping(j);end
row={n.position_mm,status,i,j,of,nf,100*(nf/of-1),oz,nz,100*(nz-oz),100*(nz/oz-1),mac,mac0,d,cost,co,cn,o.status,n.status};
end

function [R,R0]=residues(r)
fit=r.selected.fit;nm=numel(fit.frequency_hz);R=complex(zeros(16,nm));R0=R;
for m=1:nm
    p=roots([fit.a(m,2),fit.a(m,1),1]);ii=find(imag(p)>0);assert(isscalar(ii),'Need stable underdamped conjugate mode.');
    q=p(ii);wr=2*pi*fit.reference_hz(m);sp=wr*q;
    v=wr*(fit.coefficients(2*m-1,:)+q*fit.coefficients(2*m,:))/(fit.a(m,1)+2*fit.a(m,2)*q);
    R0(:,m)=v.';R(:,m)=v.'.*exp(-sp*r.delay_matrix_seconds(:));
end
assert(all(isfinite(R(:))));
end

function v=MAC(a,b)
den=sum(abs(a).^2)*sum(abs(b).^2);
if den<=realmin,v=0;else,v=min(1,max(0,abs(a'*b)^2/den));end
end

function H=transfer(r,f)
fit=r.selected.fit;f=f(:);nm=numel(fit.reference_hz);q=1i*f./fit.reference_hz(:).';
A=1+q.*fit.a(:,1).'+q.^2.*fit.a(:,2).';B=complex(zeros(numel(f),2*nm+double(fit.direct)));
for m=1:nm,B(:,2*m-1)=1./A(:,m);B(:,2*m)=q(:,m)./A(:,m);end
if fit.direct,B(:,end)=1;end
H=(B*fit.coefficients).*exp(-1i*2*pi*f.*r.delay_matrix_seconds(:).');assert(all(isfinite(H(:))));
end

function L=transform(H,Ty,Tu)
L=complex(zeros(size(H)));
for k=1:size(H,1),G=Ty*reshape(H(k,:),4,4)*Tu;L(k,:)=G(:).';end
end

function ii=bandMask(f,b)
ii=f>=b(1)&f<b(2);if b(2)==800,ii=f>=b(1)&f<=800;end
assert(nnz(ii)>=2,'Insufficient band samples.');
end

function v=complexRMS(f,H)
assert(numel(f)==size(H,1)&&numel(f)>=2&&all(diff(f)>0));
v=sqrt(trapz(f,sum(abs(H).^2,2))/(f(end)-f(1)));
end

function v=fieldOr(s,name,fallback)
if isfield(s,name),v=s.(name);else,v=fallback;end
end
