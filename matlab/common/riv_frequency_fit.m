function [iterations,beta,info]=riv_frequency_fit(f,G,weight,beta0,degrees,origins,tolerance,maxIterations,bounds)
% Frequency-domain additive RIV, paper Eq. (38), with column scaling.
% Extract only the block-diagonal coefficients after every update.
assert(isvector(f) && isreal(f) && all(isfinite(f(:)) & f(:)>0) && ...
    ismatrix(G) && ~isempty(G) && size(G,2)==numel(f) && all(isfinite(G(:))), ...
    'identification:InvalidFrequencyData', ...
    'RIV requires finite positive frequencies and matching finite FRFs.');
f=f(:); noni=size(G,1); K=numel(degrees);
sizes=2+noni*(degrees+1); offsets=[0,cumsum(sizes)]; p=sum(sizes);
if nargin<9 || isempty(bounds), bounds=[zeros(K,1),inf(K,1)]; end
assert(isequal(size(beta0),[p,K]),'RIV parameter size mismatch.');
beta=beta0; [X,Z,Y,prediction]=regression(beta);
residual=apply_sqrt_weight(stack_complex(G-prediction),weight);
cost=residual.'*residual; costs=cost; iterations=0; converged=false; reason='maximum_iterations';
for it=1:maxIterations
    Xw=apply_sqrt_weight(X,weight); Zw=apply_sqrt_weight(Z,weight);
    Yw=apply_sqrt_weight(Y,weight);
    sx=max(vecnorm(Xw),realmin); sz=max(vecnorm(Zw),realmin);
    moment=(Zw./sz).'*(Xw./sx);
    proposal=lsqminnorm(moment,(Zw./sz).'*Yw,1e-12)./sx.';
    projected=zeros(size(beta));
    for mode=1:K
        ix=offsets(mode)+(1:sizes(mode)); b=proposal(ix,mode);
        poles=roots([b(2),b(1),1]);
        poles=complex(-abs(real(poles)),imag(poles));
        den=real(poly(poles)); den=den/den(end);
        b(1:2)=[den(2);den(1)]; projected(ix,mode)=b;
    end
    accepted=false;
    for backtrack=0:20
        fraction=2^-backtrack; candidate=beta+fraction*(projected-beta);
        valid=true;
        for mode=1:K
            ix=offsets(mode)+(1:2); a=candidate(ix,mode);
            hz=1/(2*pi*sqrt(a(2)));
            valid=valid && all(isfinite(a)) && all(a>0) && hz>=bounds(mode,1) && hz<=bounds(mode,2);
        end
        if ~valid, continue; end
        [Xnew,Znew,Ynew,Pnew]=regression(candidate);
        rnew=apply_sqrt_weight(stack_complex(G-Pnew),weight);
        newCost=rnew.'*rnew;
        if newCost<=cost*(1+1e-12), accepted=true; break; end
    end
    iterations=it;
    if ~accepted, reason='no_descent'; break; end
    change=norm((candidate-beta).*sx.','fro')/max(norm(beta.*sx.','fro'),eps);
    beta=candidate; X=Xnew; Z=Znew; Y=Ynew; cost=newCost;
    costs(end+1,1)=cost; %#ok<AGROW>
    if change<tolerance
        if fraction>=0.5
            converged=true; reason='relative_parameter_change';
        else
            reason='stalled_line_search';
        end
        break;
    end
end
info=struct('converged',converged,'reason',reason,'iterations',iterations,'cost_history',costs(:));

    function [X,Z,Y,pred]=regression(parameters)
        N=numel(f); X=zeros(2*N*noni,p); Z=X; Y=zeros(2*N*noni,K); pred=zeros(noni,N);
        for frequency=1:N
            s=1i*2*pi*f(frequency); A=zeros(1,K); P=zeros(noni,K);
            for submodel=1:K
                indices=offsets(submodel)+(1:sizes(submodel)); localBeta=parameters(indices,submodel);
                A(submodel)=1+localBeta(1)*s+localBeta(2)*s^2;
                coeff=reshape(localBeta(3:end),noni,degrees(submodel)+1);
                P(:,submodel)=coeff*(s.^(0:degrees(submodel))).'/(s^origins(submodel)*A(submodel));
            end
            pred(:,frequency)=sum(P,2); x=zeros(noni,p); z=x; y=zeros(noni,K);
            for submodel=1:K
                indices=offsets(submodel)+(1:sizes(submodel)); residualPlant=G(:,frequency)-pred(:,frequency)+P(:,submodel);
                numerator=kron(s.^(0:degrees(submodel)),eye(noni))/(s^origins(submodel)*A(submodel));
                x(:,indices)=[-s*residualPlant/A(submodel),-s^2*residualPlant/A(submodel),numerator];
                z(:,indices)=[-s*P(:,submodel)/A(submodel),-s^2*P(:,submodel)/A(submodel),numerator];
                y(:,submodel)=residualPlant/A(submodel);
            end
            rows=(frequency-1)*2*noni+(1:2*noni);
            X(rows,:)=[real(x);imag(x)]; Z(rows,:)=[real(z);imag(z)]; Y(rows,:)=[real(y);imag(y)];
        end
    end
end
