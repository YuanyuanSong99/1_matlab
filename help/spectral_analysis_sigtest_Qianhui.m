%%%%%%%%%%%%%%%%%%%简单的功率谱分析+检验%%%%%%%%%%%%%%%%%%%
%导入数据（以温度为例）
  tss = siad_0';
  temp_rq_summer=tss(1,1:end-1);
%%%%%%%%%%%%%%%%%%对原始温度序列进行傅里叶变换，求频谱%%%%%%%%%%%%%%%%%%%%%%
%温度的fft
  fs=1;%fs--sample frequency:3h 1次,8次/天
  n=length(temp_rq_summer);%数据长度（采样点数）
  frequency=(fs/n)*(0:1:n-1);%频率分辨率:f_div=fs/n
%做傅里叶变换，看一下PM2.5和温度的周期
  X_temp_rqsounding=fft(temp_rq_summer);
%power,sum(ampiltude^2)计算各点的模
  X_power_temp_rqsounding=abs(X_temp_rqsounding);
%计算真实幅值：fft第一个点是直流信号，模是原始信号的N倍，其他点模是原始的N/2倍
  X_power_temp_rqsounding(1)=X_power_temp_rqsounding(1)./n;
  X_power_temp_rqsounding(2:end)=X_power_temp_rqsounding(2:end)./(n/2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%对原始温度序列计算一阶自回归%%%%%%%%%%%%%%%%%%%%%%%%%%
clear yy xx  
yy=temp_rq_summer;
  xx=tss(1,2:end);%X往后错开1个点
%用多项式拟合法，y=p1X^n+p2X^(n-1)+…+pn系数p是按照降幂排序的
  [pp]=polyfit(xx,yy,1);%RQ
%也可以用多元线性回归的方法，矩阵第一列是常数，第二列是x1，然后x2……
%X=[ones(size(xx))' xx'];p=regress(yy',X);
%生成100组随机数
for i=1:100
   %均值是0，方差是1的随机数 
   noise=randn(length(yy),1);
   %y的每一列是一组测试序列,这些新序列斜率一样，只是加了白噪声
   %就是为了检验原来频谱的peak是不是随机产生的
   yy_test(:,i)=pp(1).*xx'+noise;
end
%对y_test做功率谱
  fs_test=1;
  n_test=length(yy_test(:,1));
  frequency_test=(fs_test/n_test)*(0:1:n_test-1);%频率分辨率:f_div=fs/n
%做傅里叶变换
  X_test_temp_rq=fft(yy_test);%会对每一列进行fft
%power,sum(ampiltude^2)计算各点的模
  X_power_test_temp_rq=abs(X_test_temp_rq);
%计算真实幅值：fft第一个点是直流信号，模是原始信号的N倍，其他点模是原始的N/2倍
  X_power_test_temp_rq(1,:)=X_power_test_temp_rq(1,:)./n_test;
  X_power_test_temp_rq(2:end,:)=X_power_test_temp_rq(2:end,:)./(n_test/2);  
%画一下测试组的频谱
  [a,b]=size(X_power_test_temp_rq);
  freq=zeros(a,b);
for i=1:b
    freq(:,i)=frequency_test;
end
figure(1);clf;
  l1=plot(freq(2:n_test/2,:),X_power_test_temp_rq(2:n_test/2,:),'-','color',[0.6 0.6 0.6],'linewidth',1.5);
  hold on
  l2=plot(frequency(2:n/2),X_power_temp_rqsounding(2:n/2),'r-','linewidth',1.5);
  set(gca,'fontsize',14,'linewidth',2);
  ylabel('Amplitude');
%   set(gca,'xlim',[0 3],'xtick',0:0.5:4)%,'xticklabel',0:0.5:4);
%   set(gca,'ylim',[0 7],'ytick',0:2:10);
%   m=legend([l2],'RQ','location','northeast');set(m,'box','off');
%   xlabel('frequency（cycles/day）');
set(gcf,'position',[0 0 800 300]);
% saveas(gcf,'f:\spectrum_temp_summer_sigtest','tif');
% clear m l1 l2

