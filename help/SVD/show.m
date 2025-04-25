clc;
clear;
% 
% 
load('svd.mat');
% ---normalized and correlation coefficient
[nt,kmod] = size(A);
for i = 1:kmod
    x(:,i) = (A(:,i)-ones(nt,1)*mean(A(:,i))) / std(A(:,i),1);
    y(:,i) = (B(:,i)-ones(nt,1)*mean(B(:,i))) / std(B(:,i),1);
    % % temp = corrcoef(A(:,i),B(:,i));
    % % r(i) = temp(1,2)
    r(i) = sum(x(:,i).*y(:,i)) / nt;
end;
% ---»æÍ¼²¢±£´æ
for i = 1:kmod
  % 
  figure();
  % ---SST
  subplot(3,1,[1]);
  contourf(lon1,lat1,squeeze(U(:,:,i)));
  caxis([-0.06 0.06]);
  colorbar;
  set(gca,'XTick',[125:5:160]);
  set(gca,'XTickLabel',{'125E' '130E' '135E' '140E' '145E' '150E' '155E'});
  set(gca,'YTick',[15:5:30]);
  set(gca,'YTickLabel',{'15N' '20N' '25N'});
  title(strcat('s',num2str(i),'(SST)'));
  % ---SSH
  subplot(3,1,[2]);
  contourf(lon1,lat1,squeeze(V(:,:,i)));
  caxis([-0.06 0.06]);
  colorbar;
  set(gca,'XTick',[125:5:160]);
  set(gca,'XTickLabel',{'125E' '130E' '135E' '140E' '145E' '150E' '155E'});
  set(gca,'YTick',[15:5:30]);
  set(gca,'YTickLabel',{'15N' '20N' '25N'});
  title(strcat('s',num2str(i),'(SSH)'));
  % ---r
  subplot(3,1,[3]);
  R = num2str(r(i));
  SCF = num2str(L(i));
  t = 1:nt;
  plot(t,x(:,i),'r',t,y(:,i),'b');
  xlabel('time');
  ylabel('normalized units');
  set(gca,'XTick',[24:60:nt]);
  set(gca,'XTickLabel',{'1990' '1995' '2000' '2005'});
  set(gca,'YTick',[-4:2:4]);
  set(gca,'YTickLabel',{'-4' '-2' '0' '2' '4'});
  title(strcat('s',num2str(i),'(SST)',' and s',num2str(i),'(SSH)', ... 
      ' scf = ',SCF(1:5),' r = ',R(1:4)));
  grid on;
  print(gcf,'-djpeg',strcat('pic\',strcat(num2str(i),'svd.jpg')));
  %
end;
