    function [Cfb num_Cfb den_Cfb] = Function_PID_HB(mass,f_bw,k)

%----------------说明------------------------
%功能：设计反馈控制器
%DATE: 2022-08-13
%输入参数：mass - 质量；
%         f_bw  - 期望的开环带宽；
%         k-控制器与被控对象之间增益；
%输出参数：Cfb     - 反馈控制器传递函数；
%         num_Cfb  - 反馈控制器传递函数分母；
%         den_Cfb - 反馈控制器传递函数分子。
%---------------END--------------------------


%------------HansButtler 给出的控制器形式------------
        M = mass/k;
        a = 3;
        Kp = M*(2*pi*f_bw)^2/a;
        f_i = f_bw/(a^2);
        f_d = f_bw/a;
        f_lp = a*f_bw;
        z_lp = 0.707; %0.5~1

        C_PID = Kp*tf([1 2*pi*f_i],[1 0])*tf([1 2*pi*f_d],[2*pi*f_d]);
        C_Lowpass = tf([(2*pi*f_lp)^2],[1 2*z_lp*2*pi*f_lp (2*pi*f_lp)^2]);

        Cfb = C_PID*C_Lowpass;

    
    [num_Cfb den_Cfb] = tfdata(Cfb,'v');

end

