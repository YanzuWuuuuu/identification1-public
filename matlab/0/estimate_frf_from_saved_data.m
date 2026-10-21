function filename=estimate_frf_from_saved_data(dataDirectory)
% Regenerate LRM FRFs from saved u1..u4/y1..y4 without a Simulink rerun.
stageDirectory=fileparts(mfilename('fullpath'));
if nargin<1, dataDirectory=stageDirectory; end
previousPath=path; cleanup=onCleanup(@()path(previousPath));
addpath(fullfile(fileparts(stageDirectory),'common'));
fs=5000; channels=4;
for k=1:channels
    input=load(fullfile(dataDirectory,sprintf('u%d.mat',k)));
    output=load(fullfile(dataDirectory,sprintf('y%d.mat',k)));
    u=input.(sprintf('u%d',k)); y=output.(sprintf('y%d',k));
    assert(isvector(u)&&isvector(y)&&numel(u)==numel(y),'Invalid saved input/output data.');
    if k==1, N=numel(u); U=zeros(channels,floor(N/2)+1); Y=U; end
    assert(numel(u)==N,'All channels must have equal sample counts.');
    fu=fft(u(:))/N; fy=fft(y(:))/N;
    U(k,:)=fu(1:floor(N/2)+1).'; Y(k,:)=fy(1:floor(N/2)+1).';
end
f=(0:floor(N/2))*fs/N; responses=cell(1,channels);
for k=1:channels
    [~,~,~,responses{k}]=lrm_miso_core(Y(k,:).',U,9,1);
end
P1=responses{1}; P2=responses{2}; P3=responses{3}; P4=responses{4};
source=struct('time_data_directory',dataDirectory,'sample_rate_hz',fs, ...
    'sample_count',N,'method','scaled_QR_LRM_common_transient');
filename=fullfile(stageDirectory,'LRM_data.mat');
save(filename,'P1','P2','P3','P4','f','source');
fprintf('Stage 0: saved %d-channel LRM data with %d frequency bins.\n',channels^2,numel(f));
end
