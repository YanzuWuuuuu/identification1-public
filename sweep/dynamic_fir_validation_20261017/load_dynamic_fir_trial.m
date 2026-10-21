function trial=load_dynamic_fir_trial(source)
% Offline candidate only: not validated for insertion into hardware feedback.
% source: 'nonparam' or 'srivc'. Coefficients are frozen using training only.
% Input: v_r=w_r-Fu*wN; v_N=z^-L*wN; physical motor u=Tu*v.
% Output: z=Ty*y; calibrated N=z^-L*zN-Fy*zr; first 3 outputs unchanged.
% The NN reference branch has explicit delay 2*L*Ts. The complete corrected
% NN response also contains FIR corrections; its group delay need not be constant.
% Do not remove the plant's physical delay.
% Ty/Tu use unit-norm logical port directions, not physical Rx/Ry/Z scales.
% Connect to the physical 4x4 plant G only; do not transform G_logical again.
% Existing hardware-controller gains require the corresponding unit conversion.
if nargin<1,source='nonparam';end
assert(any(strcmp(source,{'nonparam','srivc'})),'Unknown source.');
s=load(fullfile(fileparts(mfilename('fullpath')),'dynamic_fir_coefficients.mat'));
hu=s.(['hu_' source]);hy=s.(['hy_' source]);L=double(s.(['delay_' source]));
trial=struct('source',source,'Ts',s.Ts,'Tu',s.Tu,'Ty',s.Ty, ...
    'hu',hu,'hy',hy,'added_delay_samples_per_side',L, ...
    'scope','Offline held-out-data trial; no hardware closed-loop claim');
% h(k+1,j) multiplies the sample delayed by k. tf Variable z^-1 preserves it.
Fu=tf(repmat({0},3,1),repmat({1},3,1),s.Ts,'Variable','z^-1');
Fy=tf(repmat({0},1,3),repmat({1},1,3),s.Ts,'Variable','z^-1');
for j=1:3
    Fu(j,1)=tf(hu(:,j).',1,s.Ts,'Variable','z^-1');
    Fy(1,j)=tf(hy(:,j).',1,s.Ts,'Variable','z^-1');
end
q=tf([zeros(1,L),1],1,s.Ts,'Variable','z^-1');
I3=ss(eye(3),'Ts',s.Ts);z31=ss(zeros(3,1),'Ts',s.Ts);z13=ss(zeros(1,3),'Ts',s.Ts);
trial.Fu=Fu;trial.Fy=Fy;trial.q=q;
trial.Ud=[I3,-ss(Fu);z13,ss(q)];
trial.Yd=[I3,z31;-ss(Fy),ss(q)];
trial.physical_input_map=s.Tu*trial.Ud;
trial.physical_output_map=trial.Yd*s.Ty;
end
