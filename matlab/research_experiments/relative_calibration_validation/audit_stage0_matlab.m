function audit_stage0_matlab()
% Read-only source/data check. Does not run the original Simulink script.
here=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(fileparts(here)));
src=fileread(fullfile(root,'matlab','0','Multiple_main_Control.m'));
names={'B1','C1','B2','C2','kesai','f'};
v=struct();
for k=1:numel(names)
    token=regexp(src,['(?m)^\s*' names{k} '\s*=\s*(\[[^\]]+\])'],'tokens','once');
    v.(names{k})=str2num(token{1}); %#ok<ST2NM> Local numeric matrix literals only.
end
names={'Jrx','Jry','mz'}; inertia=zeros(1,3);
for k=1:3
    token=regexp(src,['(?m)^\s*' names{k} '\s*=\s*([\d.]+)'],'tokens','once');
    inertia(k)=str2double(token{1});
end
s=tf('s');
P=v.C1*diag(1./inertia)*v.B1/s^2;
for k=1:5
    P=P+(v.C2(:,k)*v.B2(k,:))/(s^2+2*v.kesai(k)*2*pi*v.f(k)*s+(2*pi*v.f(k))^2);
end
hz=[2 20 100 v.f 800];
fr=freqresp(P,2*pi*hz);
out=struct('hz',hz,'frf_real',real(fr),'frf_imag',imag(fr));
U=zeros(20000,4);
for k=1:4
    q=load(fullfile(root,'matlab','0',sprintf('u%d.mat',k)));
    U(:,k)=q.(sprintf('u%d',k));
end
out.saved_input_channel_rms=sqrt(mean(U.^2,1));
out.saved_input_total_rms=sqrt(mean(sum(U.^2,2)));
[~,~,V]=svd(v.B1); n=V(:,end);
out.saved_input_null_energy_fraction=sum((U*n).^2)/sum(U.^2,'all');
out.saved_input_sample_count=size(U,1);
% DFT band power using the script's 5kHz sampling; DC/Nyquist excluded.
uf=fft(U)/size(U,1); freq=(0:size(U,1)-1)'*5000/size(U,1);
for lo=[2 20]
    idx=freq>=lo & freq<=800;
    out.(sprintf('saved_input_band_%d_800_rms',lo))=sqrt(2*sum(abs(uf(idx,:)).^2,'all'));
end
out.provenance='Existing stage-0 saved simulation inputs; generation history not verified, not hardware data.';
outdir=fullfile(here,'perturbation_outputs');
fid=fopen(fullfile(outdir,'audit_matlab.json'),'w');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(out));
disp(rmfield(out,{'frf_real','frf_imag'}));
end
