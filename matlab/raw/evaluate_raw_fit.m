function evaluation=evaluate_raw_fit(u,y,fit,ntrain)
% Calibrate initial conditions using TRAINING outputs only, then free simulate.
% Holdout outputs affect reported errors but never the model or prediction.
ts=fit.sample_time; no=size(y,2);
[forced,initial]=raw_modal_basis(u,ts,fit.frequency_hz,fit.damping,fit.hold_method);
prediction=forced*reshape(fit.residues,no,[]).';
ic=scaled_lsq(initial(1:ntrain,:),y(1:ntrain,:)-prediction(1:ntrain,:));
prediction=prediction+initial*ic;
train=1:ntrain; validation=ntrain+1:size(y,1);
residual=y-prediction; metrics=zeros(1,6);
metrics(1)=norm(residual(train,:),'fro')/norm(y(train,:),'fro');
metrics(2)=norm(residual(validation,:),'fro')/norm(y(validation,:),'fro');
bands={[5,200],[200,720],[720+1e-6,800],[200,800]};
for k=1:numel(bands)
    errorBand=raw_residual_transform(residual(validation,:),ts,bands{k});
    dataBand=raw_residual_transform(y(validation,:),ts,bands{k});
    metrics(k+2)=norm(errorBand,'fro')/max(norm(dataBand,'fro'),realmin);
end
evaluation=struct('prediction',prediction,'initial_coefficients',ic,'metrics',metrics, ...
    'per_output_validation_error',vecnorm(residual(validation,:))./vecnorm(y(validation,:)), ...
    'validation_samples',validation,'metric_names',{{'train_time','holdout_time', ...
    'holdout_5_200','holdout_200_720','holdout_over720_800','holdout_200_800'}});
end
