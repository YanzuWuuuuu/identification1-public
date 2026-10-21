function summary=summarize_d_residual()
%SUMMARIZE_D_RESIDUAL Aggregate the locked 4 x 2 x 2 residual-D experiment.
folder=fileparts(mfilename('fullpath')); outdir=fullfile(folder,'d_residual_outputs');
scenarios={'original','near','mid','far'}; rows=struct([]); tail=struct([]);
for ss=1:4
    for seed=1:2
        pair=cell(2,1);
        for includeD=0:1
            tag=sprintf('%s_seed%02d_D%d',scenarios{ss},seed,includeD);
            d=jsondecode(fileread(fullfile(outdir,[tag '.json'])));
            pair{includeD+1}=load(fullfile(outdir,[tag '.mat']),'Y','sigmaMean','f','core','truth','fit');
            q=pair{includeD+1};
            assert(all(q.core.D(:)==0)&&all(q.truth.D(:)==0));
            assert(numel(q.core.lambda)==5 && isequal(q.core.R,q.truth.R(:,:,1:5)));
            if ~includeD, assert(all(q.fit.estimate.D(:)==0)&&q.fit.num_free_parameters==134); end
            r.scenario=scenarios{ss}; r.seed=seed; r.include_D=includeD; r.file=[tag '.json'];
            r.cost=d.solver.mean_whitened_residual;
            r.frequency_error_percent=100*max(abs(d.frequency_relative_error));
            r.damping_error_percent=100*max(abs(d.damping_relative_error));
            r.dynamic_error_percent=100*max(d.relative.output_max_error,d.relative.input_max_error);
            r.full_band_frf_error_percent=100*d.fit.full_object_relative_frf_l2_error;
            r.flex_band_frf_error_percent=100*d.fit.band_200_800_hz.full_object_relative_frf_l2_error;
            r.rigid_to_target_five=d.baseline.joint.leak_to_flex_amplitude_ratio;
            ret=arrayfun(@(x)x.joint.retention_amplitude_ratio,d.mode);
            r.max_retention_deviation_percent=100*max(abs(ret-1));
            r.min_retention_percent=100*min(ret); r.max_retention_percent=100*max(ret);
            r.high_tail_to_target_five=d.baseline.joint.high_mode_tail_to_target_five_amplitude_ratio;
            r.estimated_D_norm=d.estimated_D_frobenius_norm;
            r.exitflag=d.solver.exitflag; r.iterations=d.solver.iterations;
            r.active_bound_count=numel(d.solver.active_bounds);
            r.normalized_gradient=d.solver.normalized_gradient;
            r.jacobian_error=d.solver.jacobian_check;
            r.free_parameters=d.solver.num_free_parameters;
            if isempty(rows), rows=r; else, rows(end+1)=r; end %#ok<AGROW>
            if seed==1 && includeD==0
                t.scenario=scenarios{ss}; t.extra_frequencies_hz=d.synthetic_extra_frequencies_hz;
                t.tail_to_five_flexible_percent=100*d.tail_diagnostic.tail_to_target_flexible_l2;
                t.D0_relative_tail_percent=100*empty_nan(d.tail_diagnostic.D0_relative_tail_l2);
                t.best_constant_relative_tail_percent=100*empty_nan(d.tail_diagnostic.best_real_constant_relative_tail_l2);
                t.upper_edge_D0_relative_tail_percent=100*empty_nan(d.tail_diagnostic.tail_relative_error_at_upper_edge);
                t.D0_norm=norm(d.tail_diagnostic.D0,'fro');
                if isempty(tail), tail=t; else, tail(end+1)=t; end %#ok<AGROW>
            end
        end
        assert(isequal(pair{1}.Y,pair{2}.Y)&&isequal(pair{1}.sigmaMean,pair{2}.sigmaMean)&&isequal(pair{1}.f,pair{2}.f));
    end
end
summary.generated=datestr(now,'yyyy-mm-dd HH:MM:SS'); summary.cases=rows; summary.tail=tail;
summary.case_count=numel(rows); summary.paired_data_verified=true;
summary.true_D_zero_verified=true; summary.original_five_residues_unchanged_verified=true;
summary.model_note='Added high modes are synthetic. True D is zero in all data. The fitted constant is only a residual approximation.';
summary.comparison='16 prespecified fits: four objects, two noise seeds, same-data D=0 vs D-free, expanded channel parameterization, 400 iterations.';
T=struct2table(rows); writetable(T,fullfile(outdir,'d_residual_cases.csv'));
fid=fopen(fullfile(outdir,'d_residual_summary.json'),'w','n','UTF-8'); assert(fid>=0);
fprintf(fid,'%s',jsonencode(summary,'PrettyPrint',true)); fclose(fid);
fig=figure('Visible','off','Color','w','Position',[100 100 1200 820]);
colors=[.16 .39 .69;.83 .35 .16];
subplot(2,2,1); bar([[tail(2:4).D0_relative_tail_percent].',[tail(2:4).best_constant_relative_tail_percent].']);
set(gca,'XTickLabel',{'Near','Mid','Far'}); ylabel('Error / true high-mode tail (%)');
title('Constant approximation of the actual tail'); legend({'Static limit D_0','Best real constant'},'Location','best'); grid on;
fields={'damping_error_percent','dynamic_error_percent','flex_band_frf_error_percent'};
titles={'Maximum damping error','Relative channel-response error','Full-object FRF error, 200-800 Hz'};
for k=1:3
    subplot(2,2,k+1); hold on;
    for inc=0:1
        for ss=1:4
            vals=T.(fields{k})(strcmp(T.scenario,scenarios{ss}) & T.include_D==inc);
            xpos=ss+(2*inc-1)*.13;
            scatter(xpos+[-.02;.02],vals,35,colors(inc+1,:),'filled','HandleVisibility','off');
            plot(xpos+[-.09,.09],[mean(vals),mean(vals)],'-','Color',colors(inc+1,:),'LineWidth',2,'HandleVisibility','off');
        end
        plot(NaN,NaN,'o-','Color',colors(inc+1,:),'DisplayName',sprintf('D %s',pick(inc,'free','fixed zero')));
    end
    set(gca,'XTick',1:4,'XTickLabel',{'Original 5','+Near 3','+Mid 3','+Far 3'},'YScale','log');
    ylabel('Relative error (%)'); title(titles{k}); legend('Location','best'); grid on;
end
sgtitle('True data contain modes only: paired five-mode fits with and without a residual D');
exportgraphics(fig,fullfile(outdir,'d_residual_comparison.png'),'Resolution',160); savefig(fig,fullfile(outdir,'d_residual_comparison.fig')); close(fig);
disp(T(:,{'scenario','seed','include_D','cost','frequency_error_percent','damping_error_percent','dynamic_error_percent','flex_band_frf_error_percent','rigid_to_target_five','exitflag'}));
fprintf('Verified %d cases, paired data/weights identical, true D always zero.\n',numel(rows));
end

function x=empty_nan(x)
if isempty(x), x=NaN; end
end
function s=pick(condition,a,b)
if condition, s=a; else, s=b; end
end
