function evaluation=evaluate_coupled_support_calibration(tag)
if nargin<1,tag='z0';end
base=fileparts(mfilename('fullpath'));out=fullfile(base,'coupled_support_20261016');
S=load(fullfile(out,['calibration_' tag '.mat']),'calib');c=S.calib;
S=load(c.source_predecessor,'calib');old=S.calib;f=c.f(:);n=numel(f);
P=c.nominal.Ty./vecnorm(c.nominal.Ty,2,2);Q=c.nominal.Tu./vecnorm(c.nominal.Tu,2,1);
H=reshape(c.H_validation.',4,4,[]);Z0=complex(zeros(size(H)));
for k=1:n,Z0(:,:,k)=[P;c.nominal.ny.']*H(:,:,k)*[Q,c.nominal.nu];end
fits={old.fits.joint,c.fits.joint};names={'previous_proportional','coupled_support'};
bands=[3 25;25 80;80 170;100 170;260 310;25 800];fitRows={};crossRows={};responses={};
for m=1:2
    fit=fits{m};ell=fit.delay_seconds;dy=ell(:,1)-ell(1,1);du=ell(1,:)-ell(1,1);
    assert(norm(ell-(ell(1,1)+dy+du),'fro')<1e-10);
    py=exp(-2i*pi*f*(max([0;dy])-dy).');pu=exp(-2i*pi*f*(max([0 du])-du));
    for delayOn=0:1
        Z=complex(zeros(size(H)));
        for k=1:n
            if delayOn,Cy=diag(py(k,:));Cu=diag(pu(k,:));else,Cy=eye(4);Cu=eye(4);end
            Z(:,:,k)=[P;fit.Ny.']*Cy*H(:,:,k)*Cu*[Q,fit.Nu];
        end
        responses{end+1}=struct('method',names{m},'delay_alignment',delayOn,'Z',Z); %#ok<AGROW>
        for b=1:size(bands,1)
            band=bands(b,:);ix=f>=band(1)-1e-6 & f<=band(2)+1e-6;w=weights(f(ix));
            for bl=1:3
                if bl==1,x=reshape(Z(4,1:3,ix),3,[]).';x0=reshape(Z0(4,1:3,ix),3,[]).';name='logical_to_N';
                elseif bl==2,x=reshape(Z(1:3,4,ix),3,[]).';x0=reshape(Z0(1:3,4,ix),3,[]).';name='N_to_logical';
                else,x=reshape(Z(4,4,ix),[],1);x0=reshape(Z0(4,4,ix),[],1);name='N_to_N';end
                rms=@(v)sqrt(sum(w.*sum(abs(v).^2,2))/sum(w));
                crossRows(end+1,:)={names{m},delayOn,band(1),band(2),name,rms(x),rms(x)/rms(x0)}; %#ok<AGROW>
            end
        end
    end
    for b=1:size(bands,1)
        band=bands(b,:);ix=f>=band(1)-1e-6 & f<=band(2)+1e-6;w=weights(f(ix));
        rms=@(v)sqrt(sum(w.*sum(abs(v).^2,2))/sum(w));
        fitRows(end+1,:)={names{m},band(1),band(2),100*rms(fit.H(ix,:)-c.H_training(ix,:))/rms(c.H_training(ix,:)), ...
            100*rms(fit.H(ix,:)-c.H_validation(ix,:))/rms(c.H_validation(ix,:))}; %#ok<AGROW>
    end
end
splitRows={};free=old.initial_free.selected.H;
for b=1:size(bands,1)
    band=bands(b,:);ix=f>=band(1)-1e-6 & f<=band(2)+1e-6;w=weights(f(ix));
    rms=@(v)sqrt(sum(w.*sum(abs(v).^2,2))/sum(w));
    fitRows(end+1,:)={'free_modal_reference',band(1),band(2),100*rms(free(ix,:)-c.H_training(ix,:))/rms(c.H_training(ix,:)), ...
        100*rms(free(ix,:)-c.H_validation(ix,:))/rms(c.H_validation(ix,:))}; %#ok<AGROW>
    splitRows(end+1,:)={band(1),band(2),100*rms(c.H_training(ix,:)-c.H_validation(ix,:))/rms(c.H_validation(ix,:))}; %#ok<AGROW>
end
splitMetrics=cell2table(splitRows,'VariableNames',{'fmin','fmax','half_record_difference_pct'});
fitMetrics=cell2table(fitRows,'VariableNames',{'method','fmin','fmax','training_relative_error_pct','validation_relative_error_pct'});
crossMetrics=cell2table(crossRows,'VariableNames',{'method','delay_alignment','fmin','fmax','block','rms','over_B0'});
angle=@(a,b)acosd(min(1,abs(a.'*b)/(norm(a)*norm(b))));
evaluation=struct('f',f,'fit_metrics',fitMetrics,'cross_metrics',crossMetrics,'split_metrics',splitMetrics,'responses',{responses}, ...
    'null_direction_change_deg',[angle(fits{1}.Ny,fits{2}.Ny),angle(fits{1}.Nu,fits{2}.Nu)], ...
    'delay_change_frobenius_us',1e6*norm(fits{2}.delay_seconds-fits{1}.delay_seconds,'fro'), ...
    'weighted_training_error_pct',[fits{1}.train_weighted_error_pct,fits{2}.train_weighted_error_pct], ...
    'rho',fits{2}.rho,'actual_support_poles',fits{2}.support_poles, ...
    'scope','One nested restricted damping-family start at z=0; same-record second-half LRM, not hardware or unbiased-model certification');
save(fullfile(out,['evaluation_' tag '.mat']),'evaluation','-v7.3');
writetable(fitMetrics,fullfile(out,['fit_comparison_' tag '.csv']));writetable(crossMetrics,fullfile(out,['cross_comparison_' tag '.csv']));
writetable(splitMetrics,fullfile(out,['split_comparison_' tag '.csv']));
disp(fitMetrics);disp(crossMetrics(crossMetrics.fmin==25 & crossMetrics.fmax==80 & ~strcmp(crossMetrics.block,'N_to_N'),:));
disp(evaluation.null_direction_change_deg);disp(evaluation.weighted_training_error_pct);
end
function w=weights(f)
df=diff(f);w=[df(1);df(1:end-1)+df(2:end);df(end)]/2;
end
