%%%%%%%%%%%%%%%%%%%%
% SAR图像频谱图获取 %
%%%%%%%%%%%%%%%%%%%%
% GBInsar 处理脚本文件
clear all;
%读取文件路径
fileName1 ='E:/GBSar/test/11/2013.07.27 - 11.39.37-StaL-000005-testnch3.dat';
fileName2 ='E:/GBSar/test/113/2013.08.02 - 11.47.20-StaL-000006-testnch3.dat';
%数据读取
row    = 2600;  
col    = 802;
spectr_x1      = zeros(1, col/2);
spectr_x2     = spectr_x1;
spectr_y1      = zeros(row, 1);
spectr_y2     = spectr_y1;
[Gsar,Phsar] = read_slcsar(fileName1,row,col);
[Hsar,Phsar] = read_slcsar(fileName2,row,col);
for i =1:row
x1=Gsar(i,:);
X=fftshift(fft(x1));
spectr_x1 =spectr_x1 +X;
x2=Hsar(i,:);
X=fftshift(fft(x2));
spectr_x2 =spectr_x2 +X;
end
for j =1:col/2;
y1=Gsar(:,j);
Y=fftshift(fft(y1));
spectr_y1 =spectr_y1 +Y;
y2=Hsar(:,j);
Y=fftshift(fft(y2));
spectr_y2 =spectr_y2 +Y;
end
labx=linspace(-150,150,col/2);
plot(labx,abs(spectr_x1));grid on;hold on;plot(labx,abs(spectr_x2),'r');
 title('距离向频谱分布');xlabel('距离向频率/MHz');
 Legend('主影像','从影像');
 labx=linspace(-4,4,row);
figure;plot(labx,abs(spectr_y1));grid on;hold on;plot(labx,abs(spectr_y2),'r');
 title('方位向频谱分布');xlabel('方位向频率/Hz');
 Legend('主影像','从影像');
