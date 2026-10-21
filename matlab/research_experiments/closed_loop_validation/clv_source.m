function [p,nominal,config]=clv_source()
% Numeric source only; never execute the stage-0 script.
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(fileparts(here)));
source=fullfile(root,'matlab','0','Multiple_main_Control.m'); txt=fileread(source);
for cellname={'B1','B2','C1','C2','f','kesai'}
    name=cellname{1}; token=regexp(txt,['(?m)^\s*' name '\s*=\s*\[([^\]]+)\]'],'tokens','once');
    rows=strsplit(token{1},';'); value=[];
    for r=1:numel(rows)
        row=sscanf(regexprep(rows{r},'[\r\n,]',' '),'%f').';
        if ~isempty(row), value=[value;row]; end %#ok<AGROW>
    end
    data.(name)=value;
end
inertia=zeros(1,3); names={'Jrx','Jry','mz'};
for j=1:3
    token=regexp(txt,['(?m)^\s*' names{j} '\s*=\s*([-+\d.eE]+)'],'tokens','once');
    inertia(j)=str2double(token{1});
end
p.M=data.C1*diag(1./inertia)*data.B1; p.D=zeros(4);
w=2*pi*data.f(:); zeta=data.kesai(:); p.lambda=-zeta.*w+1i*w.*sqrt(1-zeta.^2);
for i=1:5, p.R(:,:,i)=data.C2(:,i)*data.B2(i,:)/(p.lambda(i)-conj(p.lambda(i))); end
p.filter_y=[20;31;9;25]*1e-6; p.filter_u=[17;11;23;35]*1e-6;
p.ty=[0;40;80;15]*1e-6; p.tu=[25;65;35;100]*1e-6;
p.a=p.filter_y+p.filter_u.'; p.b=p.filter_y*p.filter_u.';
% Nominal geometric errors are fixed and only used for controller / comparator.
rng(712,'twister'); Cn=data.C1; Bn=data.B1;
for j=1:3
    Cn(:,j)=Cn(:,j)+.01*norm(Cn(:,j))*randn(4,1)/2;
    Bn(j,:)=Bn(j,:)+.01*norm(Bn(j,:))*randn(1,4)/2;
end
nominal.M=Cn*diag(1./inertia)*Bn;
config=struct('source',source,'Fs',5000,'Ts',1/5000,'period_samples',5000, ...
    'frequencies_hz',(20:4:800).','validation_frequencies_hz',(22:4:798).', ...
    'force_peak_per_line',.1,'sensor_noise_std',100e-9,'burn_samples',10000, ...
    'pilot_groups',32,'main_groups',32,'starts',16,'max_iterations',400, ...
    'max_function_evaluations',4000,'frequency_prior_hz',[200 800], ...
    'damping_bounds',[.001 .08],'filter_a_bounds',[2e-6 300e-6], ...
    'filter_b_bounds',[1e-12 (150e-6)^2], 'delay_bounds',[-100e-6 250e-6], ...
    'total_delay_bounds',[0 400e-6],'filter_region_alpha',1000, ...
    'filter_region_beta',600,'source_frequencies_hz',data.f,'D_is_fixed_zero',true);
end
