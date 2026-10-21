function realData=stack_complex(data)
% Frequency-major ordering: [real(G(:,1));imag(G(:,1));...].
realData=reshape([real(data);imag(data)],[],1);
end
