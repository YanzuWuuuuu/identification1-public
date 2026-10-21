function x=scaled_lsq(A,b)
scale=vecnorm(A);
scale=max(scale,max(max(scale)*1e-14,realmin));
x=lsqminnorm(A./scale,b,1e-12)./scale.';
end
