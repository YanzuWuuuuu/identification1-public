function results = grounding_ratio(outputDir)
%GROUNDING_RATIO Audit rigid/flexible ratios from the explicit model in 0.
% No legacy scripts, Simulink, FE solves, or control-system toolbox are run.
% Ratios use physical displacement/force ports, with equal independent input
% spectra for integrated RMS values. They are not trajectory-independent
% output ratios for arbitrary forcing. The five modes in 0 are all retained.
here = fileparts(mfilename('fullpath'));
matlabDir = fileparts(fileparts(here));
if nargin < 1, outputDir = fullfile(here,'results'); end
if ~isfolder(outputDir), mkdir(outputDir); end
sourceFile = fullfile(matlabDir,'0','Multiple_main_Control.m');
txt = fileread(sourceFile);
B1 = read_matrix(txt,'B1'); B2 = read_matrix(txt,'B2');
C1 = read_matrix(txt,'C1'); C2 = read_matrix(txt,'C2');
f0 = read_matrix(txt,'f');
inertias = [read_scalar(txt,'Jrx'),read_scalar(txt,'Jry'),read_scalar(txt,'mz')];
zeta0 = read_matrix(txt,'kesai');
M = C1*diag(1./inertias)*B1;
[UC,~,~] = svd(C1); [~,~,VB] = svd(B1);
ny = UC(:,end); nu = VB(:,end);
bands = [2 20;20 100;100 250;250 550;600 800;2 800];
bands = [bands;[f0(:).*(1-3*zeta0(:)),f0(:).*(1+3*zeta0(:))]];
selected = [1 2 10 20 100 f0 800];
f = unique([logspace(log10(2),log10(800),6001), ...
    linspace(2,800,40001),bands(:)',selected]);
s = 2i*pi*f;
N = zeros(4,4,5);
for k = 1:5, N(:,:,k) = C2(:,k)*B2(k,:); end
left = {eye(4),ny',eye(4),ny',UC(:,1:3)'};
right = {eye(4),eye(4),nu,nu,VB(:,1:3)};
labels = {'physical_4x4','output_null','input_null','both_null','orthonormal_rigid_ports'};
results = struct;
results.source = sourceFile;
results.definition = ['Pointwise ratio ||G_r(jw)||_F/||G_f(jw)||_F; ', ...
    'band RMS ratio sqrt(integral ||G_r||_F^2 df / integral ||G_f||_F^2 df). ', ...
    'Equal independent flat input spectra; no channel filtering or measurement noise.'];
results.primary = 'Five modes and damping from matlab/0; no D in this original plant.';
results.frequencies_hz = f0;
results.inertias_Jrx_Jry_mz = inertias;
results.M = M;
results.null_rigid_coefficient_numerical = abs(ny'*M*nu);
results.null_rigid_coefficient_exact = 0;
results.modal_both_null_retention = zeros(1,5);
for k = 1:5
    results.modal_both_null_retention(k) = abs(ny'*N(:,:,k)*nu)/norm(N(:,:,k),'fro');
end
results.stl = stl_properties(fullfile(matlabDir,'anlaysis','old_wafer_matlab.STL'),2700);
results.stl.ratio_to_0_Jrx_Jry_mz = ...
    [results.stl.inertia_CoM(1:2),results.stl.mass]./inertias;
pointwise = struct;
for iz = 1:2
    if iz == 1, zeta = zeta0; name = 'primary_0_zeta_001';
    else, zeta = 0.001*ones(1,5); name = 'sensitivity_zeta_0001'; end
    gr = zeros(4,4,numel(f)); gf = gr;
    for j = 1:numel(f)
        gr(:,:,j) = M/s(j)^2;
        for k = 1:5
            w0 = 2*pi*f0(k);
            gf(:,:,j) = gf(:,:,j) + N(:,:,k)/(s(j)^2+2*zeta(k)*w0*s(j)+w0^2);
        end
    end
    sc = struct; sc.damping = zeta;
    for p = 1:numel(labels)
        nr = zeros(size(f)); nf = nr;
        for j = 1:numel(f)
            nr(j) = norm(left{p}*gr(:,:,j)*right{p},'fro');
            nf(j) = norm(left{p}*gf(:,:,j)*right{p},'fro');
        end
        ratio = nr./nf;
        ps = struct;
        ps.selected_frequencies_hz = selected;
        [~,is] = ismember(selected,f);
        ps.selected_ratios = ratio(is);
        for b = 1:size(bands,1)
            take = f>=bands(b,1) & f<=bands(b,2);
            ps.bands(b).range_hz = bands(b,:);
            ps.bands(b).rms_ratio = sqrt(trapz(f(take),nr(take).^2)/trapz(f(take),nf(take).^2));
            ps.bands(b).pointwise_minmax = [min(ratio(take)),max(ratio(take))];
        end
        if p>1 && p<5
            ps.interpretation = 'The true rigid term is exactly zero; reported residual is floating-point roundoff.';
        end
        sc.(labels{p}) = ps;
        pointwise.(name).(labels{p}) = struct('rigid_norm',nr,'flexible_norm',nf,'ratio',ratio);
    end
    sc.per_entry = struct;
    for b = 1:size(bands,1)
        take = f>=bands(b,1) & f<=bands(b,2);
        ratio = abs(gr(:,:,take))./abs(gf(:,:,take));
        sc.per_entry.bands(b).range_hz = bands(b,:);
        sc.per_entry.bands(b).rms_ratio = sqrt(trapz(f(take),abs(gr(:,:,take)).^2,3)./ ...
            trapz(f(take),abs(gf(:,:,take)).^2,3));
        sc.per_entry.bands(b).pointwise_min = min(ratio,[],3);
        sc.per_entry.bands(b).pointwise_max = max(ratio,[],3);
    end
    results.(name) = sc;
end
save(fullfile(outputDir,'grounding_ratio.mat'),'results','pointwise','f');
fid = fopen(fullfile(outputDir,'grounding_ratio.json'),'w','n','UTF-8');
assert(fid>=0,'Cannot write ratio JSON.');
cleaner = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(results,PrettyPrint=true));
clear cleaner
fig = figure('Visible','off','Color','w');
q = pointwise.primary_0_zeta_001.physical_4x4;
loglog(f,q.ratio,'LineWidth',1.3); hold on; yline(1,'--'); yline(1e5,':');
xlim([1 800]); grid on; xlabel('Frequency (Hz)'); ylabel('Rigid / flexible FRF Frobenius norm');
title('Actual five-mode model in matlab/0, damping 0.01');
exportgraphics(fig,fullfile(outputDir,'grounding_ratio.png'),'Resolution',160); close(fig);
fprintf('Saved model-grounded ratio results to %s\n',outputDir);
disp(results.primary_0_zeta_001.physical_4x4);
end

function A = read_matrix(txt,name)
tok = regexp(txt,['(?m)^\s*',name,'\s*=\s*\[([^\]]+)\]'],'tokens','once');
assert(~isempty(tok),'Missing matrix %s',name);
rows = regexp(tok{1},';','split'); A = [];
for j = 1:numel(rows)
    clean = regexprep(rows{j},'[\r\n,]',' ');
    vals = sscanf(clean,'%f')';
    if ~isempty(vals), A = [A;vals]; end %#ok<AGROW>
end
end

function a = read_scalar(txt,name)
tok = regexp(txt,['(?m)^\s*',name,'\s*=\s*([-+\d.eE]+)'],'tokens','once');
assert(~isempty(tok),'Missing scalar %s',name); a = str2double(tok{1});
end

function out = stl_properties(path,rho)
txt = fileread(path);
tok = regexp(txt,'vertex\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)','tokens');
v = zeros(numel(tok),3);
for j = 1:numel(tok), v(j,:) = cellfun(@str2double,tok{j}); end
assert(mod(size(v,1),3)==0,'Invalid ASCII STL triangles.');
vol = 0; first = zeros(1,3); second = first;
for j = 1:3:size(v,1)
    tri = v(j:j+2,:); vt = dot(tri(1,:),cross(tri(2,:),tri(3,:)))/6;
    vol = vol+vt; first = first+vt*sum(tri,1)/4;
    second = second + vt/10*(sum(tri.^2,1)+tri(1,:).*tri(2,:)+ ...
        tri(1,:).*tri(3,:)+tri(2,:).*tri(3,:));
end
assert(vol>0,'Expected positive oriented closed STL volume.');
center = first/vol; mass = rho*vol;
inertia = rho*(sum(second)-second)-mass*(sum(center.^2)-center.^2);
edges = zeros(size(v,1),2); [~,~,vid] = unique(round(v,12),'rows');
for j = 1:3:size(v,1)
    t = vid(j:j+2); edges(j:j+2,:) = sort([t([1 2 3]),t([2 3 1])],2);
end
[~,~,edgeId] = unique(edges,'rows'); counts = accumarray(edgeId,1);
assert(all(counts==2),'STL surface is not closed with two faces per edge.');
out = struct('source',path,'density',rho,'volume',vol,'mass',mass, ...
    'centroid',center,'inertia_CoM',inertia,'bounds',[min(v,[],1);max(v,[],1)], ...
    'triangles',size(v,1)/3,'closed_edge_count',numel(counts));
end
