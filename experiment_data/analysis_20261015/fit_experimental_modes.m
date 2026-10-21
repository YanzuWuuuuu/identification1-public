function fit_experimental_modes
% Shared-pole local fits to preliminary MIMO spectra, not final calibration.
out=fileparts(mfilename('fullpath'));
positions=[-.5,0,.5,1];
windows={[240,315],[315,385],[455,535],[650,735]};
bounds={[260,300],[320,348;348,375],[475,510],[660,692;695,720]};
records=cell(numel(positions),1);rows={};groupRows={};
for ip=1:numel(positions)
    tag=sprintf('z%g',positions(ip));
    D=load(fullfile(out,['auto_identify_motor1234_' tag '_spectral.mat']));
    groups=cell(4,1);
    for ig=1:4
        win=windows{ig};lim=bounds{ig};
        train=localFit(D.train,win,lim,[]);
        validation=localFit(D.validation,win,lim,train.frequency_hz);
        % Sensitivity to a modest fit-window change, using the same training data.
        narrow=localFit(D.train,win+[4,-4],lim,train.frequency_hz);
        f=D.validation.f;idx=f>=win(1)&f<=win(2);
        measured=reshape(D.validation.H(:,:,idx),16,[]).';
        pred=evaluateLocal(train,f(idx));
        omegaScale=(f(idx)/mean(f(idx))).^2;
        valerr=norm((pred-measured).*omegaScale,'fro')/norm(measured.*omegaScale,'fro');
        groups{ig}=struct('train',train,'validation',validation,'window_sensitivity',narrow, ...
            'validation_relative_error',valerr,'validation_f',f(idx),'validation_data',measured,'validation_prediction',pred);
        groupRows(end+1,:)={positions(ip),ig,win(1),win(2),100*train.relative_error,100*valerr,train.exitflag}; %#ok<AGROW>
        for im=1:numel(train.frequency_hz)
            rows(end+1,:)={positions(ip),ig,im,train.frequency_hz(im), ...
                100*train.damping(im),validation.frequency_hz(im),100*validation.damping(im), ...
                abs(train.frequency_hz(im)-validation.frequency_hz(im)), ...
                abs(train.frequency_hz(im)-narrow.frequency_hz(im)), ...
                train.residue_rank1_energy(im),train.halfpower_width_hz(im), ...
                100*valerr,train.exitflag}; %#ok<AGROW>
        end
        fprintf('%s group%d f=%s damping(%%)=%s training=%.3f%% validation=%.3f%%\n',tag,ig, ...
            mat2str(train.frequency_hz.',6),mat2str(100*train.damping.',4),100*train.relative_error,100*valerr);
    end
    records{ip}=struct('position_mm',positions(ip),'groups',{groups});
end
modalTable=cell2table(rows,'VariableNames',{'z_mm','group','component','frequency_hz','damping_percent', ...
    'second_half_frequency_hz','second_half_damping_percent','split_frequency_difference_hz', ...
    'window_frequency_difference_hz','rank1_residue_energy_fraction','estimated_width_hz', ...
    'group_validation_error_percent','solver_exitflag'});
groupTable=cell2table(groupRows,'VariableNames',{'z_mm','group','fmin_hz','fmax_hz','training_error_percent','validation_error_percent','solver_exitflag'});
writetable(modalTable,fullfile(out,'local_modal_estimates.csv'));
writetable(groupTable,fullfile(out,'local_fit_validation.csv'));
save(fullfile(out,'local_modal_models.mat'),'records','modalTable','groupTable','windows','bounds');
disp(modalTable);

% Shared scale comparison across positions and between halves of each record.
ref=load(fullfile(out,'auto_identify_motor1234_z0_spectral.mat'));
bands=[25 100;240 315;315 385;455 535;650 735;25 800];comparisonRows={};
for ip=1:numel(positions)
    D=load(fullfile(out,sprintf('auto_identify_motor1234_z%g_spectral.mat',positions(ip))));
    for b=1:size(bands,1)
        id=ref.whole.f>=bands(b,1)&ref.whole.f<=bands(b,2);
        omega=(ref.whole.f(id)/mean(ref.whole.f(id))).^2;
        A=reshape(D.whole.H(:,:,id),16,[]).'.*omega;
        B=reshape(ref.whole.H(:,:,id),16,[]).'.*omega;
        T=reshape(D.train.H(:,:,id),16,[]).'.*omega;
        V=reshape(D.validation.H(:,:,id),16,[]).'.*omega;
        comparisonRows(end+1,:)={positions(ip),bands(b,1),bands(b,2), ...
            100*norm(A-B,'fro')/norm(B,'fro'),100*norm(T-V,'fro')/norm((T+V)/2,'fro')}; %#ok<AGROW>
    end
end
comparison=cell2table(comparisonRows,'VariableNames',{'z_mm','fmin_hz','fmax_hz','difference_from_z0_percent','within_record_split_difference_percent'});
writetable(comparison,fullfile(out,'position_comparison.csv'));

% Six identified local pole candidates; do not force the simulation's five modes.
fig=figure('Visible','off','Color','w','Position',[60 60 1400 900]);
tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
names={'Broad mode near 280 Hz','Additional feature near 335 Hz','Mode near 360 Hz', ...
    'Mode near 492 Hz','Additional feature near 680 Hz','Mode near 705 Hz'};
for m=1:6
    ax=nexttile(tl);hold(ax,'on');
    ix=m:6:height(modalTable);z=modalTable.z_mm(ix);ft=modalTable.frequency_hz(ix);
    fv=modalTable.second_half_frequency_hz(ix);
    plot(ax,z,ft,'o-','LineWidth',1.5,'MarkerSize',6,'DisplayName','First 30 s fit');
    plot(ax,z,fv,'s--','LineWidth',1.2,'MarkerSize',5,'DisplayName','Second 30 s fit');
    grid(ax,'on');xticks(ax,positions);xlabel(ax,'z position (mm)');ylabel(ax,'Fitted frequency (Hz)');
    title(ax,names{m});if m==1,legend(ax,'Location','best');end
end
title(tl,'Position dependence: shared-pole local fits; z = -1 mm excluded due to abnormal motor command');
exportgraphics(fig,fullfile(out,'position_modal_frequencies.png'),'Resolution',170);close(fig);

% Reference-position 4x4 Bode magnitudes and phases: curves only in fitted bands.
R=records{2};D=ref;f=D.validation.f;ib=f>=25&f<=800;
for phasePlot=[false,true]
    fig=figure('Visible','off','Color','w','Position',[40 40 1500 1000]);
    tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','compact');
    for y=1:4
        for u=1:4
            ax=nexttile(tl);hold(ax,'on');h=squeeze(D.validation.H(y,u,ib));
            if phasePlot,obs=unwrap(angle(h))*180/pi;else,obs=20*log10(max(abs(h),realmin));end
            semilogx(ax,f(ib),obs,'Color',[.65 .65 .65],'LineWidth',.7,'DisplayName','Second-half FRF');
            for ig=1:4
                g=R.groups{ig};ff=g.validation_f;hh=g.validation_prediction(:,y+4*(u-1));
                if phasePlot
                    model=unwrap(angle(hh))*180/pi;
                    % Only add an integer 360-degree branch to the model for display.
                    target=interp1(f(ib),obs,ff(1));model=model+360*round((target-model(1))/360);
                else
                    model=20*log10(max(abs(hh),realmin));
                end
                semilogx(ax,ff,model,'Color',[.05 .32 .65],'LineWidth',1.25,'HandleVisibility','off');
            end
            set(ax,'XScale','log');grid(ax,'on');xlim(ax,[25 800]);xticks(ax,[25 100 300 800]);
            title(ax,sprintf('u_%d to y_%d',u,y));
            if y==4,xlabel(ax,'Frequency (Hz)');end
            if u==1
                if phasePlot,ylabel(ax,'Phase (deg)');else,ylabel(ax,'Magnitude (dB, raw units)');end
            end
        end
    end
    if phasePlot
        title(tl,'z = 0 mm: second-half data (gray), first-half local modal fits (blue)');
        fn='z0_local_bode_phase.png';
    else
        title(tl,'z = 0 mm: held-out time segment (gray), local parametric models (blue); no fit below 25 Hz');
        fn='z0_local_bode_magnitude.png';
    end
    exportgraphics(fig,fullfile(out,fn),'Resolution',170);close(fig);
end
end

function fit=localFit(D,win,limits,start)
id=D.f>=win(1)&D.f<=win(2);f=D.f(id);
G=reshape(D.H(:,:,id),16,[]).';
nm=size(limits,1);omegaScale=(f/mean(f)).^2;
normalizer=norm(G.*omegaScale,'fro');
if isempty(start)
    indicator=sqrt(sum(abs(G).^2,2)).*f.^2;
    start=zeros(nm,1);
    for i=1:nm
        sel=find(f>=limits(i,1)&f<=limits(i,2));
        [~,a]=max(movmedian(indicator(sel),3));start(i)=f(sel(a));
    end
end
lb=log([limits(:,1);repmat(.00015,nm,1)]);ub=log([limits(:,2);repmat(.12,nm,1)]);
options=optimoptions('lsqnonlin','Display','off','MaxIterations',160, ...
    'MaxFunctionEvaluations',2500,'FunctionTolerance',1e-11,'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-8,'FiniteDifferenceType','central');
best=Inf;
for initialDamping=[.003,.015,.04]
    x0=log([start(:);repmat(initialDamping,nm,1)]);
    [xx,cost,~,flag,output]=lsqnonlin(@objective,min(max(x0,lb),ub),lb,ub,options);
    if cost<best,best=cost;x=xx;bestFlag=flag;bestOutput=output;end
end
[~,coef,pred]=objective(x);
fn=exp(x(1:nm));damp=exp(x(nm+1:end));rankFraction=zeros(nm,1);
for i=1:nm
    residue=reshape(coef(2*i-1,:)+1i*coef(2*i,:),4,4);
    sig=svd(residue);rankFraction(i)=sig(1)^2/sum(sig.^2);
end
fit=struct('frequency_hz',fn,'damping',damp,'coefficient',coef,'window_hz',win, ...
    'frequency_limits',limits,'center_hz',mean(f),'scale_hz',std(f), ...
    'frequency_hz_grid',f,'prediction',pred,'relative_error',sqrt(best), ...
    'exitflag',bestFlag,'iterations',bestOutput.iterations, ...
    'residue_rank1_energy',rankFraction,'halfpower_width_hz',2*damp.*fn);

    function [r,coef,pred]=objective(x)
        X=design(f,exp(x(1:nm)),exp(x(nm+1:end)),mean(f),std(f));
        A=[real(X.*omegaScale);imag(X.*omegaScale)];
        B=[real(G.*omegaScale);imag(G.*omegaScale)];
        cs=vecnorm(A);coef=(A./cs)\B;coef=coef./cs.';
        pred=X*coef;r=(A*coef-B)/normalizer;r=r(:);
    end
end

function G=evaluateLocal(fit,f)
G=design(f,fit.frequency_hz,fit.damping,fit.center_hz,fit.scale_hz)*fit.coefficient;
end

function X=design(f,fn,damp,center,scale)
s=2i*pi*f(:);X=complex(zeros(numel(f),2*numel(fn)+4));
for j=1:numel(fn)
    pole=2*pi*fn(j)*(-damp(j)+1i*sqrt(1-damp(j)^2));
    a=1./(s-pole);b=1./(s-conj(pole));
    X(:,2*j-1)=a+b;X(:,2*j)=1i*(a-b);
end
v=(f(:)-center)/scale;
X(:,end-3:end)=[ones(size(v)),1i*ones(size(v)),v,1i*v];
end
