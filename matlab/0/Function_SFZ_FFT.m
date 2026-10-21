function [Freq_hz output] = Function_SFZ_FFT(input,fs)

%----------------说明------------------------
%功能：对input进行FFT，返回各频率分量对应的复数
%DATE: 2016-10-30
%Author：SONG FAZHI
%Modify：SONG FAZHI 20230313
%输入参数：input-要进行傅里叶变换的序列；
%         fs-采样频率（获取input的采样频率，该频率应大于2倍最大频率分量的频率）
%输出参数：F  - FFT变换后各个分量的频率；
%         output  - FFT变换后各个分量(复数)；
%---------------END--------------------------


% 
% N = length(input);             %数据点数
% 
% y = fft((input));       %消除序列中的线性趋势项 detrend
% 
% N0 = round(N/2);                 %FFT后的结果是对称的，取一半
% Freq_hz = 0:fs/N:(N0-1)*fs/N;    %FFT变换后各个分量的频率
% y = y(1:N0);
% 
% yy(1) = abs(y(1))/N;             %直流分量(注：如果input经过了detrend，那么直流分量部分将不存在)
% yy(2:N0) = y(2:N0)/(N/2);        %各交流分量
% 
% output = yy;



    N = length(input); %数据点数
    n = N;
    Y = fft(input,n);

    if mod(n,2) == 0   %偶数
        N0 = round(n/2);
        f = fs/n*(0:N0);

        Y = Y/n;
        Y_out = Y(1:N0+1);
        Y_out(2:end-1) = 2*Y_out(2:end-1);
    else         %奇数
        N0 = round(n/2);
        f = fs/n*(0:N0-1);

        Y = Y/n;
        Y_out = Y(1:N0);
        Y_out(2:end) = 2*Y_out(2:end);
    end

    Freq_hz = f;
    output = Y_out;


end


