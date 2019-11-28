% GBInsar ����ű��ļ�
clear all;
%��ȡ�ļ�·��
path_data ='E:/GBSar/test/DataPart/';
fileSar = dir(fullfile(path_data,'*.dat'));
temp  = size(fileSar);
filenum= temp(1);


dis_los_t= zeros(1,filenum-1); %ԭʼ������λ����
dis_los2_t=zeros(1,filenum-1); %ȥ��������λ������λ����
dis_los_b= zeros(1,filenum-1); %ԭʼ������λ����
dis_los2_b=zeros(1,filenum-1); %ȥ��������λ������λ����
x1=155;
y1=255; %PSC��Ȥ��λ��
x2=385;
y2=222;
%���ݶ�ȡ
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
APS_phase  =  SumG; %������λ��
%PSC������ֵ
num = 1; %PSC����
num_h=[];%�к�
num_l=[]; %�к�
tic;
%��һ�α���������롢��λɢָ����ͳ�Ʋ���
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
Mean_G=SumG/(filenum-1);    %��ȡƽ������ͼ
Mean_coh = Sum_coh/(filenum-1);
Mean_phase = Sum_phase/(filenum-1); %��ȡƽ����λ
plot_fanmap(120,Mean_coh);      %����ƽ�����ͼ
%�ڶ��α�����������λ�Ƽ�ȥ��������λ��ֵ
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
        APS_phase = S_pha - M_pha - Mean_phase;  %������λ����ȡ
        APS_phase = avg_filter (APS_phase,3);
        APS_phase = APS_phase + Mean_phase;                     
        Sum_los2 =  Sum_los2 + ((S_pha - M_pha - APS_phase)*17.4)/(4*pi);
        Interphase = S_phsar - M_phsar - APS_phase;  %��ȡȥ������λ��ĸ�����λ
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

%ȷ��PS��
for m = 1:row
    for n= 1:col/2
       if PSC_gstd(m,n) > 0.25 || Mean_coh(m,n) < 0.8 || Mean_delta(m,n) > 0.8 %�˼�����Ե��������ֵѡȡ
        PSC_gstd(m,n) = NaN;
        PSC_mask(m,n) = NaN;
       else
           PSC_gstd(m,n) = 100 / PSC_gstd(m,n) ;
           PSC_mask(m,n) = 1;
           num_h(num) = m;                %  ��¼PSC�����к�
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
 ylabel('�ػ��״�ʵ�������봹�����ݽϲ�/mm');xlabel('�������');
 Legend('δ��������������ϲ�','��������������ϲ�');

toc;