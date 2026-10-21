function audit = audit_hundred_hz_positions
% Diagnostic only: no raw-data or fitted-model changes.
base=fileparts(mfilename('fullpath'));
parent=fullfile(base,'position_spring_identification_20261016');
out=fullfile(parent,'hundred_hz_source_audit');
if ~isfolder(out),mkdir(out);end
positions=[-.5 0 .5 1]; components=[1 2 3 4 13 5 6 8 9];
Q=[1 1 1 1;1 1 -1 -1;1 -1 -1 1;1 -1 1 -1]'/2;
assert(norm(Q'*Q-eye(4),'fro')<1e-12);
rows={}; Rall=complex(zeros(4,4,numel(components),4));
Uall=complex(zeros(4,numel(components),4));Vall=Uall;
for k=1:4
    tag=strrep(strrep(sprintf('%g',positions(k)),'-','m'),'.','p');
    S=load(fullfile(parent,['model_z' tag '.mat']));r=S.result;fit=r.selected.fit;
    assert(fit.converged && strcmp(r.status,'converged'));
    for j=1:numel(components)
        m=components(j);a=fit.a(m,:);wr=2*pi*fit.reference_hz(m);
        qp=(-a(1)+1i*sqrt(4*a(2)-a(1)^2))/(2*a(2));p=wr*qp;
        Rd=reshape((fit.coefficients(2*m-1,:)+fit.coefficients(2*m,:)*qp)/((a(1)+2*a(2)*qp)/wr),4,4);
        R=Rd.*exp(-p*r.delay_matrix_seconds);
        [U,Sv,V]=svd(R);sv=diag(Sv);u=U(:,1);v=V(:,1);
        [~,pivot]=max(abs(u));rot=exp(-1i*angle(u(pivot)));u=u*rot;v=v*rot;
        E=sum(abs(Q'*R).^2,2)/norm(R,'fro')^2;
        Ed=sum(abs(Q'*Rd).^2,2)/norm(Rd,'fro')^2;
        Rall(:,:,j,k)=R;Uall(:,j,k)=u;Vall(:,j,k)=v;
        rows(end+1,:)={positions(k),m,fit.frequency_hz(m),fit.damping(m),sv(1)^2/sum(sv.^2),sv(2)/sv(1),E(1),E(2),E(3),E(4),Ed(4),norm(R,'fro')}; %#ok<AGROW>
    end
end
metrics=cell2table(rows,'VariableNames',{'position_mm','component','frequency_hz','damping_ratio','rank1_energy_fraction','sigma2_over_sigma1','z_pattern_energy_fraction','roll_pattern_energy_fraction','pitch_pattern_energy_fraction','twist_pattern_energy_fraction','delay_removed_twist_energy_fraction','residue_norm'});
writetable(metrics,fullfile(out,'residue_pattern_metrics.csv'));
macRows={};shapeRows={};
for k=1:4
    for j=1:numel(components)
        U0=Uall(:,j,2);V0=Vall(:,j,2);r0=Rall(:,:,j,2);r1=Rall(:,:,j,k);
        umac=abs(U0'*Uall(:,j,k))^2;vmac=abs(V0'*Vall(:,j,k))^2;
        rmac=abs(r0(:)'*r1(:))^2/(norm(r0,'fro')^2*norm(r1,'fro')^2);
        cross=abs(Uall(:,[4 5 6],2)'*Uall(:,j,k)).^2;
        macRows(end+1,:)={positions(k),components(j),umac,vmac,rmac,cross(1),cross(2),cross(3)}; %#ok<AGROW>
        for ch=1:4
            shapeRows(end+1,:)={positions(k),components(j),ch,real(Uall(ch,j,k)),imag(Uall(ch,j,k)),abs(Uall(ch,j,k)),angle(Uall(ch,j,k))*180/pi,real(Vall(ch,j,k)),imag(Vall(ch,j,k))}; %#ok<AGROW>
        end
    end
end
matching=cell2table(macRows,'VariableNames',{'position_mm','component','output_MAC_to_z0_same_index','input_MAC_to_z0_same_index','residue_MAC_to_z0_same_index','output_MAC_to_z0_component4','output_MAC_to_z0_component13','output_MAC_to_z0_component5'});
shapes=cell2table(shapeRows,'VariableNames',{'position_mm','component','channel','output_real','output_imag','output_magnitude','output_phase_deg','input_real','input_imag'});
writetable(matching,fullfile(out,'position_shape_matching.csv'));
writetable(shapes,fullfile(out,'singular_vectors.csv'));
audit=struct('metrics',metrics,'matching',matching,'shapes',shapes,'observed_residues',Rall,'output_vectors',Uall,'input_vectors',Vall,'basis',Q,'positions_mm',positions,'components',components,'scope','Existing closed-loop total-command empirical fits only; no physical mode labels assigned; no new identification run. Orthogonal four-point energy fractions are geometric diagnostics, not mechanical modal energies.');
save(fullfile(out,'audit.mat'),'audit');
disp(metrics(ismember(metrics.component,[4 13 5 6]),:));
disp(matching(ismember(matching.component,[4 13 5 6]),1:5));
end
