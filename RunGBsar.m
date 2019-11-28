% GBInsar 处理脚本文件
clear all;
%读取文件路径
path_data ='E:/GBSar/test/DataPart/';
fileSar = dir(fullfile(path_data,'*.dat'));
temp  = size(fileSar);
filenum= temp(1);


dis_los_t= zeros(1,filenum-1); %原始视线向位移组
dis_los2_t=zeros(1,filenum-1); %去除大气相位视线向位移组
dis_los_b= zeros(1,filenum-1); %原始视线向位移组
dis_los2_b=zeros(1,filenum-1); %去除大气相位视线向位移组
x1=155;
y1=255; %PSC兴趣点位置
x2=385;
y2=222;
%数据读取
row    = 2600;  
col     = 802;

SumG = zeros(row, col/2);
Sum_coh = SumG;
PSC_gstd  = SumG;
PSC_mask  = SumG;
Mean_delta= SumG; 
delta      = SumG; 
Sum_los  =  SumG;
Sum_los2 = SumG;
Sum_phase  =  SumG;
APS_phase  =  SumG; %大气相位屏
%PSC点行列值
num = 1; %PSC点数
num_h=[];%行号
num_l=[]; %列号
tic;
%第一次遍历求振幅离、相位散指数等统计参数
for i = 1:filenum
    fileName = [path_data,fileSar(i).name];
    [Gsar,Phsar] = read_slcsar(fileName,row,col);
    SumG = SumG+Gsar;
    if i == 1
        M_sar = Gsar;
        M_phsar =Phsar;
    else
        S_sar = Gsar;
        S_phsar = Phsar;
        output_coh = coherence( M_sar, S_sar );
        output_coh (find (output_coh<=0.8)) = 0; 
        Sum_coh  = Sum_coh+ output_coh;
        Sum_phase =Sum_phase + S_phsar - M_phsar; 
        M_sar = S_sar;
        M_phsar = S_phsar;
    end    
end
Mean_G=SumG/(filenum-1);    %获取平均功率图
Mean_coh = Sum_coh/(filenum-1);
Mean_phase = Sum_phase/(filenum-1); %获取平均相位
plot_fanmap(120,Mean_coh);      %绘制平均相干图
%第二次遍历求视线向位移及去除大气相位均值
for i = 1:filenum-1
    fileName = [path_data,fileSar(i).name];
    [Gsar,Phsar] = read_slcsar(fileName,row,col);
    if i == 1
        M_pha = Phsar;
    else
        S_pha = Phsar;
        Sum_los = Sum_los + ((S_pha - M_pha)*17.4)/(4*pi);
        dis_los_t(i)=  Sum_los(x1,y1) ;
        dis_los_b(i)=  Sum_los(x2,y2) ;
        APS_phase = S_pha - M_pha - Mean_phase;  %大气相位屏提取
        APS_phase = avg_filter (APS_phase,3);
        APS_phase = APS_phase + Mean_phase;                     
        Sum_los2 =  Sum_los2 + ((S_pha - M_pha - APS_phase)*17.4)/(4*pi);
        Interphase = S_phsar - M_phsar - APS_phase;  %提取去大气相位后的干涉相位
        path_phase = sprintf('%s%s.mat','E:/GBSar/phasefile/',fileSar(i).name(1:42));
        save (path_phase,'Interphase');
        Mean_delta = Mean_delta + abs(S_pha - M_pha - APS_phase);
        dis_los2_t(i)= Sum_los2(x1,y1);
        dis_los2_b(i)=  Sum_los2(x2,y2) ;
        M_pha =  S_pha;
      end      
     PSC_gstd  =PSC_gstd + (Gsar - Mean_G).^2; 
end
PSC_gstd  = sqrt(PSC_gstd/(filenum-2))./Mean_G;
Mean_delta = Mean_delta/ (filenum-2); 

%确定PS点
for m = 1:row
    for n= 1:col/2
       if PSC_gstd(m,n) > 0.25 || Mean_coh(m,n) < 0.8 || Mean_delta(m,n) > 0.8 %顾及相干性的振幅离差法阈值选取
        PSC_gstd(m,n) = NaN;
        PSC_mask(m,n) = NaN;
       else
           PSC_gstd(m,n) = 100 / PSC_gstd(m,n) ;
           PSC_mask(m,n) = 1;
           num_h(num) = m;                %  记录PSC点行列号
           num_l(num)  = n;
           num = num + 1;
        end
    end
end
plot_fanmap(120,(Sum_los.*PSC_mask));
plot_fanmap(120,(Sum_los2.*PSC_mask));
save ('E:/GBSar/phasefile/mask.mat','PSC_mask');
figure;
t = 1 : filenum-1;
%plot(t,dis_los,'b',t,dis_los2,'r');
 bar(t,dis_los_b,0.4);%-dis_los_b,0.4);
 hold on;
 bar(t+0.4,dis_los2_b,0.4,'r');%-dis_los2_b,0.4,'r');
 ylabel('地基雷达实测数据与垂线数据较差/mm');xlabel('监测周期');
 Legend('未经大气改正结果较差','经大气改正结果较差');

toc;