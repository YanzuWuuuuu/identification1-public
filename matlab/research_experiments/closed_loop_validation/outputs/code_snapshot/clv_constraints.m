function [v,eq,gv,geq]=clv_constraints(x,c,cfg)
a=c.as*exp(x(c.ao+(1:16)));b=c.bs*exp(x(c.bo+(1:16)));
alpha=cfg.filter_region_alpha;beta=cfg.filter_region_beta;
dy=[0;x(c.to+(1:3))]*c.ts;du=x(c.to+(4:7))*c.ts;ell=dy+du.';
A=reshape(x(1:9),3,3);[Us,Ss,Vs]=svd(A);sigma=Ss(3,3);
v=[(2*alpha*b-a)/c.as;-(1-alpha*a+alpha^2*b); ...
    (4*b-a.^2-4*b.^2*beta^2)/(c.as^2); ...
    (cfg.total_delay_bounds(1)-ell(:))/c.ts; ...
    (ell(:)-cfg.total_delay_bounds(2))/c.ts; ...
    (c.amin-sigma)/c.amin];eq=[];
if nargout>2
    gv=zeros(c.n,numel(v));geq=[];
    for edge=1:16
        ia=c.ao+edge;ib=c.bo+edge;
        gv(ia,edge)=-a(edge)/c.as;gv(ib,edge)=2*alpha*b(edge)/c.as;
        gv(ia,16+edge)=alpha*a(edge);gv(ib,16+edge)=-alpha^2*b(edge);
        gv(ia,32+edge)=-2*a(edge)^2/c.as^2;
        gv(ib,32+edge)=(4*b(edge)-8*beta^2*b(edge)^2)/c.as^2;
        [j,k]=ind2sub([4 4],edge);ids=c.to+3+k;if j>1,ids=[ids,c.to+j-1];end
        gv(ids,48+edge)=-1;gv(ids,64+edge)=1;
    end
    gv(1:9,end)=-reshape(Us(:,3)*Vs(:,3)',9,1)/c.amin;
end
end
