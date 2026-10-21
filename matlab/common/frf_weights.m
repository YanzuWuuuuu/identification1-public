function weight=frf_weights(G,floorRelative)
% Diagonal W_ii=1/max(|G_i|,floor), NOT inverse covariance.
floorValue=max(max(abs(G),[],2)*floorRelative,eps);
w=1./max(abs(G),floorValue);
weight=reshape([w;w],[],1);
end
