function results = fea_sparse_check(outputDir)
%FEA_SPARSE_CHECK Rebuild the supplied FE geometry without dense matrices.
% Matches FEA_CB_ROM geometry, inserted vertices, material and mesh Hmin.
% Stores compact audit values only; it never edits the legacy FE model.
here = fileparts(mfilename('fullpath'));
matlabDir = fileparts(fileparts(here));
if nargin < 1, outputDir = fullfile(here,'results'); end
if ~isfolder(outputDir), mkdir(outputDir); end
stlPath = fullfile(matlabDir,'anlaysis','old_wafer_matlab.STL');
sourcePath = fullfile(matlabDir,'anlaysis','FEA_CB_ROM.m');
tAll = tic;
results = struct('source_script',sourcePath,'stl',stlPath, ...
    'status','not_started','frequency_unit','Hz','inertia_unit','kg m^2', ...
    'mass_unit','kg','origin_reference',[0 0 0]);
results.configuration = struct('youngs_modulus',70e9,'poisson_ratio',0.33, ...
    'mass_density',2700,'mesh_Hmin',0.02,'eigs_requested',20, ...
    'eigs_shift',-1,'eigs_maxit',500,'eigs_tolerance',1e-8, ...
    'flexible_frequency_cutoff_hz',1);
results.legacy_0_reference = struct('Jrx',0.0439,'Jry',0.0384,'mz',4.2809, ...
    'flexible_frequencies_hz',[286.43 344.58 488.27 705.67 738.83]);
try
    gm = fegeometry(stlPath);
    act_loadedVertex = [0 0 0; ...
        0.15*sqrt(2)/2 -0.15*sqrt(2)/2 0; ...
        -0.15*sqrt(2)/2 -0.15*sqrt(2)/2 0; ...
        0 0.15 0];
    [gm,insertedVertexIds] = addVertex(gm,Coordinates=act_loadedVertex);
    results.configuration.inserted_vertices = act_loadedVertex;
    results.configuration.inserted_vertex_ids = insertedVertexIds;
    model = femodel(AnalysisType="structuralModal",Geometry=gm);
    model.MaterialProperties = materialProperties(YoungsModulus=70e9, ...
        PoissonsRatio=0.33,MassDensity=2700);
    model = generateMesh(model,Hmin=0.02);
    FEM = assembleFEMatrices(model);
    M = FEM.M; K = FEM.K;
    assert(issparse(M) && issparse(K),'FE assembly unexpectedly returned dense matrices.');
    xyz = model.Geometry.Mesh.Nodes;
    n = size(xyz,2);
    assert(all(size(M)==3*n) && isequal(size(K),size(M)), ...
        'Unexpected structural FE degree-of-freedom ordering/dimensions.');
    results.mesh = struct('nodes',n,'dofs',size(M,1),'nnz_mass',nnz(M), ...
        'nnz_stiffness',nnz(K),'both_matrices_sparse',true, ...
        'bounds',[min(xyz,[],2)';max(xyz,[],2)']);
    qx = [ones(n,1);zeros(2*n,1)];
    qy = [zeros(n,1);ones(n,1);zeros(n,1)];
    qz = [zeros(2*n,1);ones(n,1)];
    T = [qx qy qz];
    massMatrix = T'*M*T;
    mass = qz'*M*qz;
    results.mass = struct('qz_M_qz',mass,'translation_mass_matrix',massMatrix, ...
        'translation_stiffness_matrix',T'*K*T, ...
        'ratio_to_0_mz',mass/results.legacy_0_reference.mz);
    % Each qz'*M*[0;0;coordinate] equals the corresponding first mass moment.
    firstMoment = zeros(1,3);
    for d = 1:3
        coordinateField = [zeros(2*n,1);xyz(d,:)'];
        firstMoment(d) = qz'*M*coordinateField;
    end
    center = firstMoment/mass;
    results.centroid = center;
    rOrigin = rotation_fields(xyz);
    rCenter = rotation_fields(xyz-center');
    inertiaOrigin = rOrigin'*M*rOrigin;
    inertiaCenter = rCenter'*M*rCenter;
    results.inertia = struct('about_origin',inertiaOrigin,'about_centroid',inertiaCenter, ...
        'diagonal_about_origin',diag(inertiaOrigin)', ...
        'diagonal_about_centroid',diag(inertiaCenter)', ...
        'ratio_centroid_Jx_to_0_Jrx',inertiaCenter(1,1)/results.legacy_0_reference.Jrx, ...
        'ratio_centroid_Jy_to_0_Jry',inertiaCenter(2,2)/results.legacy_0_reference.Jry);
    parallelAxis = mass*((center*center')*eye(3)-center'*center);
    results.inertia.parallel_axis_absolute_residual = ...
        norm(inertiaOrigin-inertiaCenter-parallelAxis,'fro');
    results.assembly_seconds = toc(tAll);
    results.status = 'mass_and_inertia_computed';
    % Save before the eigensolve so mass survives solver errors/interruption.
    save_compact(results,outputDir);
    tEig = tic;
    try
        opts = struct('maxit',500,'tol',1e-8,'issym',true,'isreal',true,'disp',0);
        [V,L,flag] = eigs(K,M,20,-1,opts);
        lambda = real(diag(L));
        [frequencyHz,order] = sort(sqrt(abs(lambda))/(2*pi));
        lambda = lambda(order); V = V(:,order);
        flex = find(frequencyHz>=1);
        flex = flex(1:min(14,numel(flex)));
        residual = zeros(numel(lambda),1);
        for j = 1:numel(lambda)
            kv = K*V(:,j); mv = M*V(:,j);
            residual(j) = norm(kv-lambda(j)*mv)/(norm(kv)+abs(lambda(j))*norm(mv)+eps);
        end
        results.eigensolve = struct('flag',flag,'seconds',toc(tEig), ...
            'sorted_eigenvalues',lambda','sorted_frequencies_hz',frequencyHz', ...
            'below_1_hz_count',sum(frequencyHz<1), ...
            'flexible_frequencies_hz',frequencyHz(flex)', ...
            'relative_equation_residual',residual', ...
            'flexible_relative_equation_residual',residual(flex)');
        if numel(flex)>=5
            results.eigensolve.first_five_relative_difference_to_0 = ...
                frequencyHz(flex(1:5))'./results.legacy_0_reference.flexible_frequencies_hz-1;
        end
        if flag==0, results.status = 'completed';
        else, results.status = 'eigensolve_returned_unconverged_pairs'; end
        if flag==0 && ~isempty(flex)
            try
                assert(all(lambda(flex)>0),'Flexible eigenvalues must be positive for the ROM.');
                results.fea14_mode_rom = build_z_rom(xyz,V(:,flex), ...
                    frequencyHz(flex),M,results);
            catch romError
                results.status = 'eigensolve_completed_rom_failed';
                results.rom_error = error_summary(romError);
            end
        end
    catch eigError
        results.status = 'mass_saved_eigensolve_failed';
        results.eigensolve = struct('seconds',toc(tEig),'error',error_summary(eigError));
    end
catch assemblyError
    results.status = 'assembly_or_mass_failed';
    results.error = error_summary(assemblyError);
end
results.total_seconds = toc(tAll);
save_compact(results,outputDir);
fprintf('FE sparse audit: %s; saved to %s\n',results.status,outputDir);
if isfield(results,'mass'), disp(results.mass); end
if isfield(results,'eigensolve'), disp(results.eigensolve); end
end

function rom = build_z_rom(xyz,Phi,flexHz,M,results)
% Reconstruct a separate FE-derived model at the analysis ports in test.m.
% The inserted FE vertices use another sign convention; do not substitute
% those vertices for the explicit actuator/sensor requests below.
n = size(xyz,2);
bottom = find(abs(xyz(3,:))<1e-6);
assert(~isempty(bottom),'No bottom-face nodes found.');
a = 0.15*sqrt(2)/2;
actRequested = [-a a 0;a a 0;0 -0.15 0;0 0 0];
senRequested = [-0.095 0.160 0;0.095 0.160 0; ...
    0.095 -0.160 0;-0.095 -0.160 0];
actIds = zeros(1,4); senIds = actIds;
for k = 1:4
    [~,ii] = min(vecnorm(xyz(:,bottom)-actRequested(k,:)',2,1));
    actIds(k) = bottom(ii);
    [~,ii] = min(vecnorm(xyz(:,bottom)-senRequested(k,:)',2,1));
    senIds(k) = bottom(ii);
end
actActual = xyz(:,actIds)'; senActual = xyz(:,senIds)';
assert(numel(unique(actIds))==4 && numel(unique(senIds))==4, ...
    'Nearest-node selection collapsed distinct analysis ports.');
assert(norm(imag(Phi),'fro')<=1e-10*max(norm(real(Phi),'fro'),1), ...
    'Unexpected complex FE eigenvectors for the undamped real ROM.');
Phi = real(Phi);
modeMasses = diag(Phi'*M*Phi);
assert(all(modeMasses>0),'Nonpositive modal mass.');
Phi = Phi./sqrt(modeMasses)';
modalMassGram = Phi'*M*Phi;
verticalPhi = zeros(size(Phi));
verticalPhi(2*n+(1:n),:) = Phi(2*n+(1:n),:);
verticalMassFraction = diag(verticalPhi'*M*verticalPhi)';
Bu = Phi(actIds+2*n,:)';
Cy = Phi(senIds+2*n,:);
center = results.centroid;
actRel = actActual-center; senRel = senActual-center;
% Coordinate order is [theta_x, theta_y, z_com]; omega cross r yields
% z-displacement = y*theta_x - x*theta_y + z_com.
Br = [actRel(:,2)';-actRel(:,1)';ones(1,4)];
Cr = [senRel(:,2),-senRel(:,1),ones(4,1)];
Jinv = results.inertia.about_centroid\eye(3);
mass = results.mass.qz_M_qz;
H = blkdiag(Jinv(1:2,1:2),1/mass);
rigidM = Cr*H*Br;
[Uc,~,~] = svd(Cr); [~,~,Vb] = svd(Br);
ny = Uc(:,end); nu = Vb(:,end);
nf = numel(flexHz); zeta = 0.01*ones(1,nf);
numerators = zeros(4,4,nf); residueNorm = zeros(1,nf);
projectionAmplitude = zeros(1,nf); projectionFraction = zeros(1,nf);
for k = 1:nf
    numerators(:,:,k) = Cy(:,k)*Bu(k,:);
    residueNorm(k) = norm(numerators(:,:,k),'fro');
    projectionAmplitude(k) = abs(ny'*numerators(:,:,k)*nu);
    if residueNorm(k)>0
        projectionFraction(k) = projectionAmplitude(k)/residueNorm(k);
    else
        projectionFraction(k) = NaN;
    end
end
selected = [1 2 10 20 100 flexHz(:)'];
bands = [2 20;20 100;100 250;250 550;600 800;2 800];
f = unique([logspace(log10(2),log10(800),6001), ...
    linspace(2,800,40001),selected,bands(:)']);
nr = zeros(size(f)); nflex = nr;
for j = 1:numel(f)
    s = 2i*pi*f(j); gr = rigidM/s^2; gf = zeros(4,4);
    for k = 1:nf
        w = 2*pi*flexHz(k);
        gf = gf+numerators(:,:,k)/(s^2+2*zeta(k)*w*s+w^2);
    end
    nr(j) = norm(gr,'fro'); nflex(j) = norm(gf,'fro');
end
ratio = nr./nflex;
[~,selectedIds] = ismember(selected,f);
rom = struct;
rom.identity = 'FEA14modeROM: current STL/material/mesh and test.m analysis ports; distinct from stage 0';
rom.damping_assumption = ['Uniform damping ratio 0.01 is an assumed value borrowed from stage 0, ', ...
    'not a damping estimate or an undamped FE result.'];
rom.model_equation = 'P(s)=rigid_M/s^2 + sum_k modal_numerator(:,:,k)/(s^2+2*zeta_k*w_k*s+w_k^2); D=0';
rom.mode_count = nf;
rom.frequencies_hz = flexHz(:)'; rom.damping = zeta;
rom.ports = struct('coordinate_source','matlab/anlaysis/test.m', ...
    'selection_rule','Nearest FE node restricted to abs(z)<1e-6, without interpolation', ...
    'actuator_requested_xyz',actRequested,'sensor_requested_xyz',senRequested, ...
    'actuator_node_ids',actIds,'sensor_node_ids',senIds, ...
    'actuator_actual_xyz',actActual,'sensor_actual_xyz',senActual, ...
    'actuator_distance_to_requested',vecnorm(actActual-actRequested,2,2)', ...
    'sensor_distance_to_requested',vecnorm(senActual-senRequested,2,2)');
rom.rigid_M = rigidM; rom.rigid_Br = Br; rom.rigid_Cr = Cr;
rom.rigid_inverse_mass_inertia = H;
rom.rigid_M_singular_values = svd(rigidM)';
rom.rigid_M_rank = rank(rigidM);
rom.output_null = ny'; rom.input_null = nu';
rom.rigid_null_residual = abs(ny'*rigidM*nu);
rom.flexible_B = Bu; rom.flexible_C = Cy;
rom.modal_numerators = numerators;
rom.modal_numerator_frobenius_norm = residueNorm;
rom.modal_both_null_amplitude = projectionAmplitude;
rom.modal_both_null_projection_fraction = projectionFraction;
rom.modal_mass_orthogonality_error_fro = norm(modalMassGram-eye(nf),'fro');
rom.modal_mass_gram = modalMassGram;
rom.modal_vertical_mass_fraction = verticalMassFraction;
rom.nearly_in_plane_modes = find(verticalMassFraction<1e-8);
rom.selected_frequencies_hz = selected;
rom.selected_rigid_frobenius = nr(selectedIds);
rom.selected_flexible_frobenius = nflex(selectedIds);
rom.selected_rigid_over_flexible_ratio = ratio(selectedIds);
rom.ratio_definition = ['Physical 4x4 FRF Frobenius norm ratio; integrated RMS assumes ', ...
    'equal independent flat input spectra. All retained flexible modes enter every frequency.'];
for b = 1:size(bands,1)
    take = f>=bands(b,1) & f<=bands(b,2);
    rom.bands(b).range_hz = bands(b,:);
    rom.bands(b).rms_rigid_over_flexible_ratio = ...
        sqrt(trapz(f(take),nr(take).^2)/trapz(f(take),nflex(take).^2));
    rom.bands(b).pointwise_ratio_minmax = [min(ratio(take)),max(ratio(take))];
end
rom.band_2_to_800_rms_ratio = rom.bands(6).rms_rigid_over_flexible_ratio;
rom.limitations = ['Finite mesh, nearest-node point ports, 14-mode truncation and assumed damping; ', ...
    'individual residues in a near-repeated modal cluster depend on the modal basis. ', ...
    'This ROM is a separate FE configuration and does not establish equality to stage 0.'];
end

function R = rotation_fields(xyz)
% Infinitesimal displacement omega cross position, axis order x,y,z.
n = size(xyz,2); x = xyz(1,:)'; y = xyz(2,:)'; z = xyz(3,:)';
R = [zeros(n,1),z,-y; -z,zeros(n,1),x; y,-x,zeros(n,1)];
end

function out = error_summary(err)
out = struct('identifier',err.identifier,'message',err.message);
if ~isempty(err.stack), out.first_stack_frame = err.stack(1); end
end

function save_compact(results,outputDir)
save(fullfile(outputDir,'fea_sparse_check.mat'),'results');
fid = fopen(fullfile(outputDir,'fea_sparse_check.json'),'w','n','UTF-8');
assert(fid>=0,'Cannot write FE sparse audit JSON.');
guard = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(results,PrettyPrint=true));
end
