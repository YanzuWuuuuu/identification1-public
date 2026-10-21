function audit=diagnose_antiresonance_weights
% Factorial diagnostic: weight exponent and clipping, same model and starts.
base=fileparts(mfilename('fullpath'));out=fullfile(base,'weight_diagnosis_20261016');
if ~isfolder(out),mkdir(out);end
S=load(fullfile(base,'frf_preprocessing_comparison.mat'),'comparison');C=S.comparison;
f=C.f;fn=[15 125 137 281 335 359 493 678 685 705];
zeta=[.15 .01 .01 .02 .02 .02 .02 .005 .005 .005];
bounds=[5 24;115 132;132 145;260 300;320 348;348 375;475 510;665 682;682 695;695 725];
opts=struct('fixed_modes',1,'max_iterations',150,'score_tolerance',1e-6,'verbose',false);
sources={'welch','lrm'};floors=[1e-3 1e-4 0];
models=struct('tag',{},'source',{},'exponent',{},'floor_relative',{},'fit',{},'H',{});
rows={};localRows={};floorRows={};
for e=1:2
    if e==1,H=C.welchTrain;else,H=C.lrmTrain;end
    for rho=floors
        floorValue=rho*max(abs(H),[],1);den=max(abs(H),max(floorValue,realmin));
        clipped=abs(H)<floorValue;
        for ch=1:16
            ii=f>=180 & f<=260;localAbs=abs(H(ii,ch));
            floorRows(end+1,:)={sources{e},rho,ch,mod(ch-1,4)+1,ceil(ch/4), ...
                floorValue(ch),min(localAbs),100*mean(clipped(:,ch)),100*mean(clipped(ii,ch)), ...
                floorValue(ch)/min(localAbs)}; %#ok<AGROW>
        end
        for power=1:2
            tag=sprintf('%s_p%d_floor%g',sources{e},power,rho);W=den.^(-power);W=W/median(W(:));
            fit=frf_srivc_fit(f,H,fn,zeta,bounds,W,opts);hm=fit.prediction;
            models(end+1)=struct('tag',tag,'source',sources{e},'exponent',power, ...
                'floor_relative',rho,'fit',fit,'H',hm); %#ok<AGROW>
            b=f>=200 & f<=735;
            relativeFloor=max(1e-3*max(abs(C.lrmTrain),[],1),eps);
            relScale=max(abs(C.lrmValidation),relativeFloor);
            rows(end+1,:)={tag,sources{e},power,rho,fit.converged,fit.reason,fit.iterations,fit.score, ...
                100*norm(hm-C.welchValidation,'fro')/norm(C.welchValidation,'fro'), ...
                100*norm(hm(b,:)-C.welchValidation(b,:),'fro')/norm(C.welchValidation(b,:),'fro'), ...
                100*norm(hm-C.lrmValidation,'fro')/norm(C.lrmValidation,'fro'), ...
                100*sqrt(mean(abs((hm-C.lrmValidation)./relScale).^2,'all')), ...
                max(fit.history(2:end,5)),max(fit.linear_design_condition)}; %#ok<AGROW>
            for ch=[1 6]
                ii=f>=180 & f<=260;fb=f(ii);h=hm(ii,ch);wl=C.welchValidation(ii,ch);lr=C.lrmValidation(ii,ch);
                [mv,mi]=min(abs(h));[wv,wi]=min(abs(wl));[lv,li]=min(abs(lr));
                localRows(end+1,:)={tag,ch,100*norm(h-wl)/norm(wl),100*norm(h-lr)/norm(lr), ...
                    fb(mi),20*log10(mv),fb(wi),20*log10(wv),fb(li),20*log10(lv), ...
                    20*log10(abs(h(li))/lv),angle(h(li)/lr(li))*180/pi}; %#ok<AGROW>
            end
            fprintf('%s: converged %d, score %.3g, iter %d, FRFval %.3f%%\n', ...
                tag,fit.converged,fit.score,fit.iterations,rows{end,9});
        end
    end
end
summary=cell2table(rows,'VariableNames',{'tag','source','exponent','floor_relative','converged','reason', ...
    'iterations','score','welch_error_25_800_percent','welch_error_200_735_percent', ...
    'lrm_error_25_800_percent','common_relative_lrm_error_percent','max_step_ridge','max_linear_design_condition'});
local=cell2table(localRows,'VariableNames',{'tag','channel','welch_notch_error_percent','lrm_notch_error_percent', ...
    'model_min_hz','model_min_db','welch_min_hz','welch_min_db','lrm_min_hz','lrm_min_db', ...
    'model_excess_db_at_lrm_min','phase_error_deg_at_lrm_min'});
clipping=cell2table(floorRows,'VariableNames',{'source','floor_relative','channel','output','input','absolute_floor', ...
    'min_training_magnitude_in_notch_band','clipped_all_percent','clipped_notch_percent','floor_to_notch_min_ratio'});
audit=struct('f',f,'models',models,'summary',summary,'local',local,'clipping',clipping, ...
    'welchTrain',C.welchTrain,'welchValidation',C.welchValidation, ...
    'lrmTrain',C.lrmTrain,'lrmValidation',C.lrmValidation, ...
    'weights','Quadratic W=max(abs(trainingFRF),rho*channelMax)^(-exponent); globally normalized', ...
    'selection','Factorial audit only; no configuration selected by validation error');
writetable(summary,fullfile(out,'weight_sweep_summary.csv'));
writetable(local,fullfile(out,'weight_sweep_notches.csv'));
writetable(clipping,fullfile(out,'weight_clipping_audit.csv'));
save(fullfile(out,'weight_sweep.mat'),'audit','-v7.3');
disp(summary);disp(local);
disp(clipping(ismember(clipping.channel,[1 6]) & clipping.floor_relative==1e-3,:));
end
