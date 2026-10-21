function metrics=frf_metrics(f,G,P,floorRelative)
bands={[5,720],[200,720],[720,800],[5,800]};
names={'identified_band','flexible_band','unmodeled_band','full_band'};
for k=1:numel(bands)
    ix=f>=bands{k}(1)&f<=bands{k}(2);
    if k==3, ix=f>bands{k}(1)&f<=bands{k}(2); end
    d=G(:,ix); residual=d-P(:,ix);
    w=frf_weights(d,floorRelative);
    metrics.(names{k})=struct('frequency_range_hz',bands{k}, ...
        'relative_complex_error',norm(residual,'fro')/max(norm(d,'fro'),eps), ...
        'relative_weighted_error',norm(apply_sqrt_weight(stack_complex(residual),w))/ ...
            max(norm(apply_sqrt_weight(stack_complex(d),w)),eps));
end
end
