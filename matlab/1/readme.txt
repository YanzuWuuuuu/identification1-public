其中para_emstimation为main函数
no = 4;
ni = 3;
isubsystem = 9;
% 选择所选频响范围
f_select = 5:1:1000;
% f_select = 5:5:320;
% frf导入
f_pri = [98.6,197.556,298.135,401.268,506.841,614.62,722.98,830,930];
load_data;
[frf_select,w_select] = select_frf(vecG,f_select,f,count);
这些参数需要修改

将生成的belta_init.m文件复制到第二步的文件中