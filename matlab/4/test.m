f_rho = zeros(no+ni+isubsystem*(no*ni+2),1);
phai_rbml = rho0{1}(1:no,1);
phai_rbmr = rho0{1}(no+1:no+ni,1);
M = phai_rbml*phai_rbmr';
vecM = M(:);
f_rho(1:no*ni,1) = vecM;
for i = 1:isubsystem
   real_lamda = rho0{k+1}(1,1);
   imag_lamda = rho0{k+1}(2,1);
   phai_flexl = rho0{k+1}(3:no+2,1);
   phai_flexr = rho0{k+1}(no+3:2+no+ni,1);
   lamda = real_lamda+1i*imag_lamda;
   a0_temp = lamda*lamda';
   a1_temp = -(lamda+lamda');
   a2 = 1/a0_temp;
   a1 = a1_temp/a0_temp;
   B = phai_flexl*phai_flexr';
   vecB = B(:)/a0_temp;
   B1{i} = [a1,a2,vecB']';
   f_rho(no*ni+(i-1)*(no*ni+2)+1:no*ni+i*(no*ni+2),1) = B1{i};
end

eee = f_rho-belta_est_total;