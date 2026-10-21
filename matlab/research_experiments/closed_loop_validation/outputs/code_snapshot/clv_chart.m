function [x,c,lb,ub]=clv_chart(p,cfg)
% Physical truth is never an argument in fitting; p is pilot initialization.
K=numel(p.lambda); c.K=K; c.ts=1e-5; c.as=5e-5; c.bs=1e-9;
c.fo=15; c.zo=15+K; c.ro=15+2*K; c.to=15+16*K; c.ao=c.to+7; c.bo=c.ao+16;
c.n=c.bo+16; c.ms=norm(p.M,2); [c.U,~,c.V]=svd(p.M);
x=zeros(c.n,1); core=c.U'*p.M*c.V/c.ms;
x(1:9)=reshape(core(1:3,1:3),[],1);x(10:12)=core(1:3,4);x(13:15)=core(4,1:3).';
x(c.fo+(1:K))=log(abs(p.lambda)/(2*pi));
x(c.zo+(1:K))=log(-real(p.lambda)./abs(p.lambda));
for i=1:K
    R=p.R(:,:,i); [~,pivot]=max(vecnorm(R,2,2));
    c.order(:,i)=[pivot,setdiff(1:4,pivot,'stable')].';
    right=R(pivot,:); c.rs(i)=norm(right); [~,col]=max(abs(right));
    left=R(:,col)/right(col); a=left(c.order(2:4,i)); b=right.'/c.rs(i);
    ix=c.ro+14*(i-1)+(1:14); x(ix)=[real(a);imag(a);real(b);imag(b)];
end
x(c.to+(1:3))=(p.ty(2:4)-p.ty(1))/c.ts;
x(c.to+(4:7))=(p.tu+p.ty(1))/c.ts;
x(c.ao+(1:16))=log(p.a(:)/c.as);x(c.bo+(1:16))=log(p.b(:)/c.bs);
lb=-20*ones(c.n,1);ub=20*ones(c.n,1);
lb(1:15)=x(1:15)-.3;ub(1:15)=x(1:15)+.3;
f=exp(x(c.fo+(1:K)));lo=max(200,.9*f);hi=min(800,1.1*f);
mid=(f(1:end-1)+f(2:end))/2;lo(2:end)=max(lo(2:end),mid);hi(1:end-1)=min(hi(1:end-1),mid);
lb(c.fo+(1:K))=log(lo);ub(c.fo+(1:K))=log(hi);
lb(c.zo+(1:K))=log(cfg.damping_bounds(1));ub(c.zo+(1:K))=log(cfg.damping_bounds(2));
lb(c.to+(1:7))=cfg.delay_bounds(1)/c.ts;ub(c.to+(1:7))=cfg.delay_bounds(2)/c.ts;
lb(c.ao+(1:16))=log(cfg.filter_a_bounds(1)/c.as);ub(c.ao+(1:16))=log(cfg.filter_a_bounds(2)/c.as);
lb(c.bo+(1:16))=log(cfg.filter_b_bounds(1)/c.bs);ub(c.bo+(1:16))=log(cfg.filter_b_bounds(2)/c.bs);
c.initial_frequency_hz=f;c.lb=lb;c.ub=ub;
end
