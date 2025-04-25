%%%%%%%%%%%%%%%%%%%�򵥵Ĺ����׷���+����%%%%%%%%%%%%%%%%%%%
%�������ݣ����¶�Ϊ����
  tss = siad_0';
  temp_rq_summer=tss(1,1:end-1);
%%%%%%%%%%%%%%%%%%��ԭʼ�¶����н��и���Ҷ�任����Ƶ��%%%%%%%%%%%%%%%%%%%%%%
%�¶ȵ�fft
  fs=1;%fs--sample frequency:3h 1��,8��/��
  n=length(temp_rq_summer);%���ݳ��ȣ�����������
  frequency=(fs/n)*(0:1:n-1);%Ƶ�ʷֱ���:f_div=fs/n
%������Ҷ�任����һ��PM2.5���¶ȵ�����
  X_temp_rqsounding=fft(temp_rq_summer);
%power,sum(ampiltude^2)��������ģ
  X_power_temp_rqsounding=abs(X_temp_rqsounding);
%������ʵ��ֵ��fft��һ������ֱ���źţ�ģ��ԭʼ�źŵ�N����������ģ��ԭʼ��N/2��
  X_power_temp_rqsounding(1)=X_power_temp_rqsounding(1)./n;
  X_power_temp_rqsounding(2:end)=X_power_temp_rqsounding(2:end)./(n/2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%��ԭʼ�¶����м���һ���Իع�%%%%%%%%%%%%%%%%%%%%%%%%%%
clear yy xx  
yy=temp_rq_summer;
  xx=tss(1,2:end);%X�����1����
%�ö���ʽ��Ϸ���y=p1X^n+p2X^(n-1)+��+pnϵ��p�ǰ��ս��������
  [pp]=polyfit(xx,yy,1);%RQ
%Ҳ�����ö�Ԫ���Իع�ķ����������һ���ǳ������ڶ�����x1��Ȼ��x2����
%X=[ones(size(xx))' xx'];p=regress(yy',X);
%����100�������
for i=1:100
   %��ֵ��0��������1������� 
   noise=randn(length(yy),1);
   %y��ÿһ����һ���������,��Щ������б��һ����ֻ�Ǽ��˰�����
   %����Ϊ�˼���ԭ��Ƶ�׵�peak�ǲ������������
   yy_test(:,i)=pp(1).*xx'+noise;
end
%��y_test��������
  fs_test=1;
  n_test=length(yy_test(:,1));
  frequency_test=(fs_test/n_test)*(0:1:n_test-1);%Ƶ�ʷֱ���:f_div=fs/n
%������Ҷ�任
  X_test_temp_rq=fft(yy_test);%���ÿһ�н���fft
%power,sum(ampiltude^2)��������ģ
  X_power_test_temp_rq=abs(X_test_temp_rq);
%������ʵ��ֵ��fft��һ������ֱ���źţ�ģ��ԭʼ�źŵ�N����������ģ��ԭʼ��N/2��
  X_power_test_temp_rq(1,:)=X_power_test_temp_rq(1,:)./n_test;
  X_power_test_temp_rq(2:end,:)=X_power_test_temp_rq(2:end,:)./(n_test/2);  
%��һ�²������Ƶ��
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
%   xlabel('frequency��cycles/day��');
set(gcf,'position',[0 0 800 300]);
% saveas(gcf,'f:\spectrum_temp_summer_sigtest','tif');
% clear m l1 l2

