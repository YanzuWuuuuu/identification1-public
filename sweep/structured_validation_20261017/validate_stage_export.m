function checks=validate_stage_export
% Compare the standalone export with source predictions, then test native
% Simulink delay/SS/gain wiring in an unsaved temporary model.
out=fileparts(mfilename('fullpath'));root=fileparts(fileparts(out));
target=fullfile(root,'stage_identified_z0.m');txt=fileread(target);
sourceBase=load(fullfile(root,'sweep','srivc_retest_20261016','model_z0.mat'));
sourceStructure=load(fullfile(out,'trial_z0.mat'));f=sourceBase.result.f;
rows={};
for choice={'srivc','structured'}
    for mapping={'geometric','band_static','structured'}
        current=strrep(txt,"model_choice = 'srivc';",['model_choice = ''' choice{1} ''';']);
        current=strrep(current,"decoupling_choice = 'band_static';",['decoupling_choice = ''' mapping{1} ''';']);
        runOutput=evalc(current); %#ok<NASGU>
        response=freqresp(G,2*pi*f);exported=reshape(response,16,[]).';
        if strcmp(choice{1},'srivc'),reference=sourceBase.result.selected.H;
        else,reference=sourceStructure.result.delay.H;end
        sourceError=norm(exported-reference,'fro')/norm(reference,'fro');
        assert(sourceError<1e-9);
        core=freqresp(G0,2*pi*f);aligned=freqresp(G_virtual_aligned,2*pi*f);
        expected=zeros(size(aligned));
        for k=1:numel(f),expected(:,:,k)=Ty*core(:,:,k)*Tu*exp(-2i*pi*f(k)*tau_common_aligned);end
        alignedError=norm(aligned(:)-expected(:))/norm(expected(:));assert(alignedError<1e-9);
        rawInputResponse=freqresp(G_recorded_input,2*pi*f);
        scaleError=norm(rawInputResponse(:)*input_record_multiplier-response(:))/norm(response(:));
        assert(scaleError<1e-10);
        rows(end+1,:)={choice{1},mapping{1},sourceError,alignedError,scaleError,cond(Tu),cond(Ty),max(real(eig(A)))}; %#ok<AGROW>
    end
end
checks.frequency=cell2table(rows,'VariableNames',{'model','mapping','source_relative_error','delay_alignment_relative_error', ...
    'input_scaling_relative_error','Tu_condition','Ty_condition','max_real_pole'});
% Restore and test the actual delivered defaults, within this function only.
runOutput=evalc(txt); %#ok<NASGU>
checkModel='stage_export_check_20261017';
assert(~bdIsLoaded(checkModel),'Temporary verification model name already in use.');
new_system(checkModel);cleanup=onCleanup(@()bdclose(checkModel));
mws=get_param(checkModel,'ModelWorkspace');
assignin(mws,'plant_A',A);assignin(mws,'plant_B',B);assignin(mws,'plant_C',C);assignin(mws,'plant_D',D);
assignin(mws,'plant_x0',x0);assignin(mws,'plant_Tu',Tu);assignin(mws,'plant_Ty',Ty);
assignin(mws,'plant_tau_input',tau_input);assignin(mws,'plant_tau_output',tau_output);
t=(0:1e-5:.16).';ramp=.5*(1-cos(pi*min(max((t-.01)/.02,0),1)));
uVirtual=.001*ramp.*sin(2*pi*t*[37 53 71 89]);
assignin(mws,'plant_test_input',[t,uVirtual]);
set_param(checkModel,'Solver','ode45','StopTime',num2str(t(end)), ...
    'MaxStep','1e-5','RelTol','1e-7','AbsTol','1e-10','ReturnWorkspaceOutputs','on');
add_block('simulink/Sources/From Workspace',[checkModel '/Input'], ...
    'VariableName','plant_test_input','Interpolate','on');
add_block('simulink/Math Operations/Gain',[checkModel '/Tu'], ...
    'Gain','plant_Tu','Multiplication','Matrix(K*u)');
add_block('simulink/Signal Routing/Demux',[checkModel '/SplitInput'],'Outputs','4');
add_block('simulink/Signal Routing/Mux',[checkModel '/MergeInput'],'Inputs','4');
add_block('simulink/Continuous/State-Space',[checkModel '/Core'], ...
    'A','plant_A','B','plant_B','C','plant_C','D','plant_D','X0','plant_x0');
add_block('simulink/Signal Routing/Demux',[checkModel '/SplitOutput'],'Outputs','4');
add_block('simulink/Signal Routing/Mux',[checkModel '/MergeOutput'],'Inputs','4');
add_block('simulink/Math Operations/Gain',[checkModel '/Ty'], ...
    'Gain','plant_Ty','Multiplication','Matrix(K*u)');
add_block('simulink/Sinks/To Workspace',[checkModel '/RecordedOutput'], ...
    'VariableName','plant_y_virtual','SaveFormat','Timeseries');
add_line(checkModel,'Input/1','Tu/1');add_line(checkModel,'Tu/1','SplitInput/1');
add_line(checkModel,'MergeInput/1','Core/1');add_line(checkModel,'Core/1','SplitOutput/1');
add_line(checkModel,'MergeOutput/1','Ty/1');add_line(checkModel,'Ty/1','RecordedOutput/1');
for k=1:4
    inName=['InputDelay' num2str(k)];outName=['OutputDelay' num2str(k)];
    add_block('simulink/Continuous/Transport Delay',[checkModel '/' inName], ...
        'DelayTime',sprintf('plant_tau_input(%d)',k),'InitialOutput','0');
    add_block('simulink/Continuous/Transport Delay',[checkModel '/' outName], ...
        'DelayTime',sprintf('plant_tau_output(%d)',k),'InitialOutput','0');
    add_line(checkModel,['SplitInput/' num2str(k)],[inName '/1']);
    add_line(checkModel,[inName '/1'],['MergeInput/' num2str(k)]);
    add_line(checkModel,['SplitOutput/' num2str(k)],[outName '/1']);
    add_line(checkModel,[outName '/1'],['MergeOutput/' num2str(k)]);
end
simulation=sim(checkModel);sig=simulation.plant_y_virtual;
uPhysical=uVirtual*Tu.';
yPhysical=lsim(G,uPhysical,t,x0,'foh');yExpected=yPhysical*Ty.';
yMatched=interp1(t,yExpected,sig.Time,'linear');yActual=sig.Data;
timeError=norm(yActual(:)-yMatched(:))/norm(yMatched(:));
checks.simulink=struct('model',checkModel,'duration_s',t(end),'max_step_s',1e-5, ...
    'relative_error_vs_lsim',timeError,'sample_count',numel(sig.Time), ...
    'scope','Native Tu + 4 scalar input Transport Delays + real State Space + 4 scalar output Transport Delays + Ty; unsaved temporary model');
assert(all(isfinite(yActual(:))) && timeError<2e-3,'Simulink export comparison failed: %.3g',timeError);
checks.saved_at=char(datetime('now'));checks.deliverable=target;
save(fullfile(out,'stage_export_checks.mat'),'checks');
writetable(checks.frequency,fullfile(out,'stage_export_frequency_checks.csv'));
fprintf('Six export configurations source/FRF match max %.3g\n',max(checks.frequency.source_relative_error));
fprintf('Simulink wiring vs MATLAB lsim relative error %.3g over %.2f s\n',timeError,t(end));
end
