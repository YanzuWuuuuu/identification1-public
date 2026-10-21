function [Y,J,G]=clv_mean(x,c,f,Ts,Kf,D)
nf=numel(f);wantJ=nargout>1;
if wantJ,[G,Jg]=clv_sampled(x,c,f,Ts);else,G=clv_sampled(x,c,f,Ts);end
Gp=reshape(G,4,4,nf);Id=repmat(eye(4),1,1,nf);
Sy=pagemldivide(Id+pagemtimes(Gp,Kf),Id);
SuD=pagemldivide(Id+pagemtimes(Kf,Gp),D);
Y=reshape(pagemtimes(pagemtimes(Sy,Gp),D),16,nf);
if wantJ
    J=reshape(pagemtimes(pagemtimes(reshape(Sy,4,4,nf,1), ...
        reshape(Jg,4,4,nf,c.n)),reshape(SuD,4,4,nf,1)),16*nf,c.n);
end
end
