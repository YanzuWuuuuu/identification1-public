function summary=summarize_validation()
%SUMMARIZE_VALIDATION Read saved cases without changing fits or evaluators.
% Primary cohort: common/repeated/distinct, seed 1:8, noise .003,
% N=8 independent records, one start, extra output pole = 0, both methods.
% Slices and mismatched cases are retained separately. A mismatched case's
% task_pass is undefined; core_task_pass never certifies the actual object.
folder=fileparts(mfilename('fullpath'));
outdir=fullfile(folder,'outputs'); resultdir=fullfile(folder,'results');
if ~isfolder(resultdir), mkdir(resultdir); end
files=dir(fullfile(outdir,'*.json'));
rows=struct([]); ignored={}; parseErrors={};
pattern=['^(?<scenario>[a-z_]+)_seed(?<seed>\d+)_noise(?<noise>[-+0-9.eE]+)' ...
    '_N(?<N>\d+)_(?<method>expanded|physical)_s(?<starts>\d+)' ...
    '_extra(?<extra>[-+0-9.eE]+)\.json$'];
for i=1:numel(files)
    tag=regexp(files(i).name,pattern,'names');
    if isempty(tag), ignored{end+1,1}=files(i).name; continue; end %#ok<AGROW>
    try
        d=jsondecode(fileread(fullfile(outdir,files(i).name)));
        row=case_row(d,tag,files(i).name);
        if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
    catch err
        parseErrors{end+1,1}=sprintf('%s: %s',files(i).name,err.message); %#ok<AGROW>
    end
end
assert(~isempty(rows),'No readable validation_case JSON files were found.');
T=struct2table(rows);
T=sortrows(T,{'scenario','seed','method','noise_level','records','starts','extra_pole_seconds'});
scenarios={'common','repeated','distinct'}; methods={'expanded','physical'};
primary=false(height(T),1); missing={}; duplicates={};
for s=1:numel(scenarios)
    for seed=1:8
        for m=1:numel(methods)
            mask=T.metadata_consistent & strcmp(T.scenario,scenarios{s}) & ...
                T.seed==seed & strcmp(T.method,methods{m}) & ...
                abs(T.noise_level-.003)<1e-14 & T.records==8 & ...
                T.starts==1 & T.extra_pole_seconds==0;
            key=sprintf('%s/seed%02d/%s/noise0.003/N8/s1/extra0',scenarios{s},seed,methods{m});
            if sum(mask)==1
                primary(mask)=true;
            elseif ~any(mask)
                missing{end+1,1}=key; %#ok<AGROW>
            else
                duplicates{end+1,1}=key; %#ok<AGROW>
            end
        end
    end
end
T.primary_cohort=primary;
P=T(primary,:);
groups=struct([]);
for s=1:numel(scenarios)
    for m=1:numel(methods)
        B=P(strcmp(P.scenario,scenarios{s}) & strcmp(P.method,methods{m}),:);
        g=group_summary(B,scenarios{s},methods{m});
        if isempty(groups), groups=g; else, groups(end+1)=g; end %#ok<AGROW>
    end
end
paired=struct([]); missingPairs={};
for s=1:numel(scenarios)
    for seed=1:8
        a=P(strcmp(P.scenario,scenarios{s}) & P.seed==seed & strcmp(P.method,'expanded'),:);
        b=P(strcmp(P.scenario,scenarios{s}) & P.seed==seed & strcmp(P.method,'physical'),:);
        if height(a)~=1 || height(b)~=1
            missingPairs{end+1,1}=sprintf('%s/seed%02d',scenarios{s},seed); %#ok<AGROW>
            continue;
        end
        q.scenario=scenarios{s}; q.seed=seed;
        q.dynamic_error_expanded_percent=a.relative_dynamic_max_percent;
        q.dynamic_error_physical_percent=b.relative_dynamic_max_percent;
        q.dynamic_expanded_minus_physical_percentage_points= ...
            a.relative_dynamic_max_percent-b.relative_dynamic_max_percent;
        q.leak_flex_expanded=a.actual_joint_leak_to_flex;
        q.leak_flex_physical=b.actual_joint_leak_to_flex;
        q.leak_flex_expanded_over_physical_db= ...
            20*log10(max(a.actual_joint_leak_to_flex,realmin)/max(b.actual_joint_leak_to_flex,realmin));
        q.retention_deviation_expanded_percent=a.retention_max_deviation_percent;
        q.retention_deviation_physical_percent=b.retention_max_deviation_percent;
        q.mean_cost_expanded_minus_physical=a.mean_whitened_residual-b.mean_whitened_residual;
        q.fit_seconds_expanded_over_physical=a.fit_seconds/max(b.fit_seconds,realmin);
        q.task_pass_expanded=logical(a.actual_task_pass);
        q.task_pass_physical=logical(b.actual_task_pass);
        q.exitflag_expanded=a.exitflag; q.exitflag_physical=b.exitflag;
        q.active_bound_count_expanded=a.active_bound_count;
        q.active_bound_count_physical=b.active_bound_count;
        if isempty(paired), paired=q; else, paired(end+1)=q; end %#ok<AGROW>
    end
end
summary.generated_local_time=datestr(now,'yyyy-mm-dd HH:MM:SS');
summary.primary_contract=struct('scenarios',{scenarios},'methods',{methods}, ...
    'seeds',1:8,'noise_level',.003,'independent_records',8,'starts',1, ...
    'extra_output_pole_seconds',0,'expected_case_count',48,'expected_pair_count',24);
summary.readable_case_count=height(T);
summary.primary_case_count=height(P); summary.primary_pair_count=numel(paired);
summary.primary_complete=height(P)==48 && isempty(missing) && isempty(duplicates);
summary.missing_primary_cases=missing; summary.duplicate_primary_cases=duplicates;
summary.missing_primary_pairs=missingPairs;
summary.metadata_inconsistent_files=T.file(~T.metadata_consistent);
summary.task_pass_inconsistent_files=T.file(~T.task_pass_agrees_with_recalculation);
summary.ignored_non_case_json=ignored; summary.parse_errors=parseErrors;
summary.group_statistics=groups; summary.paired_cases=paired;
summary.primary_cases=table2struct(P);
summary.other_cases=table2struct(T(~primary,:));
summary.interpretation={ ...
    'Task pass is a predeclared engineering threshold, separate from solver termination.', ...
    'Positive exit flags are solver-reported termination, not a proof of global optimality.', ...
    'A fixed 200-iteration budget is paired; elapsed time and parameter counts are not equal.', ...
    'Mismatch cases have no actual task-pass label: stored pass/dynamics/retention refer to the core.', ...
    'The actual mismatch leakage ratio is used only where separately evaluated.', ...
    'Modal retention removes the common dynamic factor; it is not total output amplitude.', ...
    'The known-error leakage bound is deterministic and is not a confidence certificate.', ...
    'Eight records and eight seeds do not establish asymptotic rates or confidence coverage.', ...
    'Whitened residuals use estimated variances; finite-N expected cost is not exactly one.'};
writetable(T,fullfile(resultdir,'validation_all_cases.csv'));
writetable(P,fullfile(resultdir,'validation_primary_cases.csv'));
writetable(struct2table(groups),fullfile(resultdir,'validation_primary_groups.csv'));
if ~isempty(paired), writetable(struct2table(paired),fullfile(resultdir,'validation_primary_pairs.csv')); end
save_json(fullfile(resultdir,'validation_summary.json'),summary);
if ~isempty(P), plot_primary(P,groups,resultdir); end
fprintf('Validation summary: %d/48 unique primary cases; %d/24 matched pairs; complete=%d.\n', ...
    height(P),numel(paired),summary.primary_complete);
fprintf('Retained separately: %d other cases; parse errors=%d; duplicate keys=%d.\n', ...
    height(T)-height(P),numel(parseErrors),numel(duplicates));
end

function r=case_row(d,t,name)
r.file=name; r.scenario=d.scenario; r.method=d.method; r.seed=d.seed;
r.noise_level=d.noise_level; r.records=d.independent_records;
r.starts=str2double(t.starts); r.extra_pole_seconds=d.extra_output_pole_seconds;
r.metadata_consistent=strcmp(d.scenario,t.scenario) && strcmp(d.method,t.method) && ...
    d.seed==str2double(t.seed) && d.independent_records==str2double(t.N) && ...
    abs(d.noise_level-str2double(t.noise))<1e-14 && ...
    abs(d.extra_output_pole_seconds-str2double(t.extra))<1e-14 && ...
    numel(d.candidate_costs)==r.starts;
r.relative_dynamic_max_percent=100*max(d.relative.output_max_error,d.relative.input_max_error);
r.frequency_max_relative_percent=100*max(abs(d.frequency_relative_error));
r.damping_max_relative_percent=100*max(abs(d.damping_relative_error));
r.rigid_matrix_relative_error=d.rigid.relative_matrix_error;
r.output_null_angle_deg=d.rigid.output_null_angle_deg;
r.input_null_angle_deg=d.rigid.input_null_angle_deg;
r.relative_delay_max_error_us=max(d.delay.output_max_abs_error_us,d.delay.input_max_abs_error_us);
r.core_joint_leak_to_flex=d.baseline.joint.leak_to_flex_amplitude_ratio;
r.core_static_leak_to_flex=d.baseline.static_only.leak_to_flex_amplitude_ratio;
r.core_nominal_leak_to_flex=d.baseline.nominal.leak_to_flex_amplitude_ratio;
r.core_joint_rigid_reduction_vs_static_db=20*log10( ...
    max(d.baseline.static_only.rigid_l2,realmin)/max(d.baseline.joint.rigid_l2,realmin));
r.core_joint_rigid_reduction_vs_nominal_db=d.baseline.joint.rigid_reduction_vs_nominal_db;
ret=arrayfun(@(x)x.joint.retention_amplitude_ratio,d.mode);
r.retention_min_ratio=min(ret); r.retention_max_ratio=max(ret);
r.retention_max_deviation_percent=100*max(abs(ret-1));
r.validation_unweighted_frf_relative_error=d.validation_relative_frf_error;
r.separability_residual=d.denominator.max_relative_separability_residual;
r.output_compensator_max_gain=d.compensation.output_max_gain;
r.input_compensator_max_gain=d.compensation.input_max_gain;
r.mean_whitened_residual=d.solver.mean_whitened_residual;
r.iterations=d.solver.iterations; r.exitflag=d.solver.exitflag;
r.solver_positive_exit=d.solver.exitflag>0;
r.active_bound_count=numel(d.solver.active_bounds);
r.first_order_optimality=d.solver.first_order_optimality;
r.normalized_gradient=d.solver.normalized_gradient;
r.jacobian_check_relative_error=d.solver.jacobian_check;
r.chart_rcond=d.solver.chart_rcond;
r.fit_seconds=d.solver.elapsed_seconds;
r.initialization_seconds=d.initialization_seconds;
r.core_task_pass=scalar_or_nan(d.task_pass);
if isfield(d,'core_task_pass'), r.core_task_pass=double(d.core_task_pass); end
r.task_pass_agrees_with_recalculation=r.core_task_pass== ...
    (r.frequency_max_relative_percent<1 && r.relative_dynamic_max_percent<5 && ...
    r.core_joint_leak_to_flex<.01 && r.retention_max_deviation_percent<10);
r.bound_gap_valid=logical(d.bound.valid_gap_condition);
r.core_oracle_bound_max_excess=scalar_or_nan(d.bound.max_absolute_excess);
r.actual_joint_leak_to_flex=r.core_joint_leak_to_flex;
r.actual_task_pass=double(r.core_task_pass);
r.actual_relative_dynamic_max_percent=r.relative_dynamic_max_percent;
r.actual_retention_max_deviation_percent=r.retention_max_deviation_percent;
r.evaluation_scope='correct_model';
if r.extra_pole_seconds>0
    r.evaluation_scope='core_metrics_except_separate_actual_leakage';
    r.actual_joint_leak_to_flex=NaN;
    if isfield(d,'mismatch') && isfield(d.mismatch,'actual_rigid_to_flexible_l2')
        r.actual_joint_leak_to_flex=d.mismatch.actual_rigid_to_flexible_l2;
    end
    r.actual_task_pass=NaN;
    r.actual_relative_dynamic_max_percent=NaN;
    r.actual_retention_max_deviation_percent=NaN;
elseif r.extra_pole_seconds<0
    r.evaluation_scope='five_true_modes_four_fitted_target_first_four';
    r.actual_task_pass=NaN;
end
end

function v=scalar_or_nan(v)
if isempty(v), v=NaN; end
end

function g=group_summary(B,scenario,method)
g.scenario=scenario; g.method=method; g.n=height(B);
g.task_pass_count=sum(B.actual_task_pass==1);
g.task_pass_rate=g.task_pass_count/max(height(B),1);
g.positive_exit_count=sum(B.solver_positive_exit);
g.zero_exit_count=sum(B.exitflag==0);
g.active_bound_case_count=sum(B.active_bound_count>0);
g.iterations_at_least_200_count=sum(B.iterations>=200);
g.dynamic_error_median_percent=med(B.relative_dynamic_max_percent);
g.dynamic_error_max_percent=mx(B.relative_dynamic_max_percent);
g.frequency_error_median_percent=med(B.frequency_max_relative_percent);
g.frequency_error_max_percent=mx(B.frequency_max_relative_percent);
g.damping_error_median_percent=med(B.damping_max_relative_percent);
g.damping_error_max_percent=mx(B.damping_max_relative_percent);
g.actual_leak_flex_median=med(B.actual_joint_leak_to_flex);
g.actual_leak_flex_max=mx(B.actual_joint_leak_to_flex);
g.retention_min_ratio=mn(B.retention_min_ratio);
g.retention_max_ratio=mx(B.retention_max_ratio);
g.retention_deviation_max_percent=mx(B.retention_max_deviation_percent);
g.rigid_reduction_vs_static_median_db=med(B.core_joint_rigid_reduction_vs_static_db);
g.rigid_reduction_vs_static_min_db=mn(B.core_joint_rigid_reduction_vs_static_db);
g.delay_error_max_us=mx(B.relative_delay_max_error_us);
g.validation_frf_error_max=mx(B.validation_unweighted_frf_relative_error);
g.normalized_gradient_median=med(B.normalized_gradient);
g.normalized_gradient_max=mx(B.normalized_gradient);
g.mean_whitened_residual_median=med(B.mean_whitened_residual);
g.fit_seconds_median=med(B.fit_seconds);
g.fit_seconds_max=mx(B.fit_seconds);
g.oracle_bound_excess_max=mx(B.core_oracle_bound_max_excess);
end

function x=med(x)
x=x(isfinite(x)); if isempty(x), x=NaN; else, x=median(x); end
end
function x=mn(x)
x=x(isfinite(x)); if isempty(x), x=NaN; else, x=min(x); end
end
function x=mx(x)
x=x(isfinite(x)); if isempty(x), x=NaN; else, x=max(x); end
end

function save_json(path,value)
fid=fopen(path,'w','n','UTF-8'); assert(fid>=0,'Could not write summary JSON.');
clean=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
end

function plot_primary(P,groups,resultdir)
fig=figure('Visible','off','Color','w','Position',[100 100 1450 800]);
clean=onCleanup(@() close(fig)); %#ok<NASGU>
colors=[.12 .38 .72;.84 .34 .12];
fields={'relative_dynamic_max_percent','actual_joint_leak_to_flex', ...
    'retention_max_deviation_percent','frequency_max_relative_percent','fit_seconds'};
titles={'Composite relative-dynamic error','Actual rigid / target-flexible L2', ...
    'Modal retention amplitude deviation','Maximum modal-frequency error','Fit elapsed time'};
ylabs={'Maximum error (%)','Amplitude ratio','Maximum over four modes (%)','Relative error (%)','Seconds'};
labels={'Common E','Common P','Repeated E','Repeated P','Distinct E','Distinct P'};
for panel=1:5
    ax=subplot(2,3,panel); hold(ax,'on');
    for g=1:numel(groups)
        B=P(strcmp(P.scenario,groups(g).scenario) & strcmp(P.method,groups(g).method),:);
        vals=B.(fields{panel}); n=numel(vals);
        if n==0, continue; end
        x=g+linspace(-.15,.15,n).'; c=colors(1+strcmp(groups(g).method,'physical'),:);
        scatter(ax,x,vals,28,c,'filled','MarkerFaceAlpha',.75);
        plot(ax,[g-.22,g+.22],[median(vals),median(vals)],'-','Color',c,'LineWidth',2);
    end
    if any(panel==[1 2 4]), set(ax,'YScale','log'); end
    set(ax,'XLim',[.5 6.5],'XTick',1:6,'XTickLabel',labels,'XTickLabelRotation',30);
    title(ax,titles{panel}); ylabel(ax,ylabs{panel}); grid(ax,'on');
end
ax=subplot(2,3,6);
counts=[[groups.positive_exit_count].',[groups.zero_exit_count].',[groups.active_bound_case_count].'];
bar(ax,counts,'grouped');
set(ax,'XTick',1:6,'XTickLabel',labels,'XTickLabelRotation',30);
title(ax,'Solver termination and bound counts'); ylabel(ax,'Cases (bound counts overlap)');
legend(ax,{'Positive exit','Exit flag 0','Any active bound'},'Location','best'); grid(ax,'on');
sgtitle(sprintf('Paired fixed-budget validation: %d/48 cases, E=expanded, P=physical; dots=seeds, bars=medians',height(P)));
exportgraphics(fig,fullfile(resultdir,'validation_primary_summary.png'),'Resolution',180);
savefig(fig,fullfile(resultdir,'validation_primary_summary.fig'));
end
