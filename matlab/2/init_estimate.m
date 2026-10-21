function [num,belta] = init_estimate(f,frf,den_para,Weight,belta0,noni,e,Max_iter,K,D_switch)

if D_switch == 0
belta1 = zeros(noni*(K+1),1);
belta2 = belta0 ;
w = 2*pi*f;
N = length(w);
G = frf;
count = 1;
phai = [];
X = [];
Y = [];
belta_count = {};
while norm(belta2-belta1)/norm(belta1)>e && count <= Max_iter
     belta1 = belta2;
     for ifrf = 1:N
         s = 1i*w(ifrf);
         R = G(:,ifrf);
         for isub = 1:K
         A = 1 + den_para{isub}(2)*s + den_para{isub}(1)*s^2;
         phai_part{isub} = (eye(noni,noni)/A)';
         phai = cat(2,phai,phai_part{isub}');
         end
         phai = [eye(noni,noni)/(s^2),phai]';
         X_part = real_ima_divide(phai');
         Y_part = real_ima_divide(R);
         X = cat(1,X,X_part);
         Y = cat(1,Y,Y_part);
         phai = [];
     end
     % belta2 = lsqminnorm(X'*Weight*X,X'*Weight*Y);
     belta2 = (X'*Weight*X)\(X'*Weight*Y);
     belta_count{count} = belta2;
     norm1(count) = norm(belta2-belta1)/norm(belta1);
     count = count +1;
     X = [];
     Y = [];
end
belta = belta_count{count-1};
num = count-1;

else
belta1 = zeros(noni*(K+2),1);
belta2 = belta0 ;
w = 2*pi*f;
N = length(w);
G = frf;
count = 1;
phai = [];
X = [];
Y = [];
belta_count = {};
while norm(belta2-belta1)/norm(belta1)>e && count <= Max_iter
     belta1 = belta2;
     for ifrf = 1:N
         s = 1i*w(ifrf);
         R = G(:,ifrf);
         for isub = 1:K
         A = 1 + den_para{isub}(2)*s + den_para{isub}(1)*s^2;
         phai_part{isub} = (eye(noni,noni)/A)';
         phai = cat(2,phai,phai_part{isub}');
         end
         phai = [eye(noni,noni)/(s^2),phai,eye(noni,noni)]';
         X_part = real_ima_divide(phai');
         Y_part = real_ima_divide(R);
         X = cat(1,X,X_part);
         Y = cat(1,Y,Y_part);
         phai = [];
     end
     % belta2 = lsqminnorm(X'*Weight*X,X'*Weight*Y);
     belta2 = (X'*Weight*X)\(X'*Weight*Y);
     belta_count{count} = belta2;
     norm1(count) = norm(belta2-belta1)/norm(belta1);
     count = count +1;
     X = [];
     Y = [];
end
belta = belta_count{count-1};
num = count-1;
end

