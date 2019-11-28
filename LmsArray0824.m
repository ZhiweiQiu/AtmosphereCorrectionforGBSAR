function [WpointR,BpointR]=LmsArray0824(Baseimg,Warpimg,Bpoint,Wpoint,Num_it,r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LMS Ϊ��С����Ӱ��ƥ�亯��
%%%%% ���������BaseimgΪ����Ӱ��WarpimgΪ��У��Ӱ��Bpoint,Wpoint�ֱ�Ϊ��׼Ӱ��ʹ�ƥ��Ӱ���������,[N X 2]
%%%%% Num_itΪ����������
%%%%% rΪƥ�䴰�ڰ뾶��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% WpointR\BpointRΪƥ�������





%Num_it=10;
%r=7;

if (size(Bpoint,1)==0)
    error('�����������������Ϊ�� ');
end

if (size(Bpoint,1)~=size(Wpoint,1))
    error('������� ');
end

BpointR=zeros(size(Bpoint,1),2);
WpointR=zeros(size(Wpoint,1),2);
BW=zeros(size(Bpoint,1),1);


if (size(size(Baseimg),2)>2)
    Baseimg=rgb2gray(Baseimg);
end
if (size(size(Warpimg),2)>2)
    Warpimg=rgb2gray(Warpimg);
end

[Bimagesizer,Bimagesizec]=size(Baseimg);
[Wimagesizer,Wimagesizec]=size(Warpimg);

Baseimg=double(Baseimg);
Warpimg=double(Warpimg);
%%%%%%%%%%%%%%%%%%��˹�˲�
BA=fspecial('gaussian',[11,11]);
Baseimg=filter2(BA,Baseimg);
Warpimg=filter2(BA,Warpimg);

for num_Bpoint=1:size(Bpoint,1)
    
    
    Lo=zeros(Num_it+1,1);
    %Lo(1)=0;
    
%     if (size(size(Baseimg),2)>2)
%         Baseimg=rgb2gray(Baseimg);
%     end
%     if (size(size(Warpimg),2)>2)
%         Warpimg=rgb2gray(Warpimg);
%     end
%     
%     Baseimg=double(Baseimg);
%     Warpimg=double(Warpimg);
%     %%%%%%%%%%%%%%%%%%��˹�˲�
%     BA=fspecial('gaussian',[11,11]);
%     Baseimg=filter2(BA,Baseimg);
%     Warpimg=filter2(BA,Warpimg);
    
    
    %%%%%%%%%%%%%%%%%%��С����ƥ��ǰ���
    Wpr=Wpoint(num_Bpoint,1);
    Wpc=Wpoint(num_Bpoint,2);
    Bpr=Bpoint(num_Bpoint,1);
    Bpc=Bpoint(num_Bpoint,2);
    
       
    if(Wpr-r-1<0||Wpr+r+1>Wimagesizer||Wpc-r-1<0||Wpc+r+1>Wimagesizec)
       %error('�������Wpc���� ');
       BW(num_Bpoint)=1;
       continue;
    end
    if(Bpr-r-1<0||Bpr+r+1>Bimagesizer||Bpc-r-1<0||Bpc+r+1>Bimagesizec)%%%%%%%%%%%%%%%%�ĳɵ�������ֱ�߾���0825
       %error('�������Wpc���� ');
       BW(num_Bpoint)=1;
       continue;
    end
    
    %%%%%%%%%%%%%%%%���ó�ֵ
    h0=0;
    h1=1;
    a0=Bpr-Wpr;
    a1=1;
    a2=0;
    b0=Bpc-Wpc;
    b1=0;
    b2=1;
    
    %%%%%%%%%%%%%%%%�������
    hh0=zeros(Num_it+1,1);
    hh0(1)=h0;
    hh1=zeros(Num_it+1,1);
    hh1(1)=h1;
    aa0=zeros(Num_it+1,1);
    aa0(1)=a0;
    aa1=zeros(Num_it+1,1);
    aa1(1)=a1;
    aa2=zeros(Num_it+1,1);
    aa2(1)=a2;
    bb0=zeros(Num_it+1,1);
    bb0(1)=b0;
    bb1=zeros(Num_it+1,1);
    bb1(1)=b1;
    bb2=zeros(Num_it+1,1);
    bb2(1)=b2;
    
    %%%%%%%%%%%%%%%%%%%��ȡ��ʼ�飬��ƥ���Ϊ���ģ���rΪ���ڰ뾶
    for i=1:2*r+3
        for j=1:2*r+3
            wx(i)=Wpr-r-2+i;
            wy(j)=Wpc-r-2+j;
            Warp1(i,j)=bilinc(Warpimg,wx(i),wy(j));
            %             Bx=Bpr-r-2+i;
            %             By=Bpc-r-2+i;
            %             BaseBxBy(i,j)=bilinc(Baseimg,Bx,By);
        end
    end
    
    for i=2:2*r+2
        for j=2:2*r+2
            WarpBxBy(i-1,j-1)=bilinc(Warpimg,wx(i),wy(j));
        end
    end
    
    
    
    %     BaseBxBy1=BaseBxBy;
    
    %     Lo(1)=Loxiangguan(Warp1,BaseBxBy);%%%%%%%%%%%�����ʼ��֮�����ϵ��
    Lo(1)=0.6;
    
    %%%%%%%%%%%%%%%%%%%%%�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:Num_it
        
        for i=1:2*r+3
            for j=1:2*r+3
                Bx=aa0(k)+aa1(k)*wx(i)+aa2(k)*wy(j);
                By=bb0(k)+bb1(k)*wx(i)+bb2(k)*wy(j);
%                 if( Bx<1||Bx>size(Baseimg,1)||Bx==NaN||By<1||By>size(Baseimg,2)||By==NaN)
%                     break;
%                 end
                BaseBxBy1(i,j)=bilinc(Baseimg,Bx,By);
            end
        end
        
        
        
        
        for i=2:2*r+2
            for j=2:2*r+2
                BaseBxBy(i-1,j-1)=hh0(k)+hh1(k)*BaseBxBy1(i,j);
            end
        end
        
        
        Lo(k+1)=Loxiangguan(WarpBxBy,BaseBxBy);%%%%%%%%%%%%%%%%%%%%%%%%������β������ĺ����ϵ��
        
        if (Lo(k+1)<=Lo(k)&&Lo(k)>=0.95)
            break;
        end
        if (Lo(k+1)<=Lo(k)&&Lo(k)<0.95)%%%%%%%%%%%%%%%%%%%�޸�һ�£��ж�Lo//0825
            %BW(num_Bpoint)=1;
            break;
        end
        c=zeros((2*r+1)^2,8);%%%%%%%%%%%%%%%%%%%%%%���ó�ʼC�����С��������ο������Ӱ������С����Ӱ��ƥ��̲ĵ�
        L=zeros((2*r+1)^2,1);%%%%%%%%%%%%%%%%%%%%%%���ó�ʼL�����С��l=g1-g2
        t=1;
        for i=2:2*r+2
            for j=2:2*r+2
                gy=(BaseBxBy1(i,j+1)-BaseBxBy1(i,j-1))/2;
                gx=(BaseBxBy1(i+1,j)-BaseBxBy1(i-1,j))/2;
                g2=BaseBxBy1(i,j);
                g1=Warp1(i,j);
                %c(t,:)=[1;g2;hh1(k)*gx;hh1(k)*i*gx;hh1(k)*j*gx;hh1(k)*gy;hh1(k)*i*gy;hh1(k)*j*gy]';
                c(t,:)=[1;g2;gx;wx(i)*gx;wy(j)*gx;gy;wx(i)*gy;wy(j)*gy]';
                L(t)=g1-g2;
                t=t+1;
            end
        end
        
        LX=inv(c'*c)*(c'*L);%%%%%%%%%%%%%%%%%%%%%����X=[dh0,dh1,da0,da1,da2,db0,db1,db2]
        
        %%%%%%%%%%%%%%%%%%%%%%%������β���%%%%%%%%%%%%%%%%%%%%%%
        hh0(k+1)=hh0(k)+LX(1)+hh0(k)*LX(2);
        hh1(k+1)=hh1(k)+hh1(k)*LX(2);
        aa0(k+1)=aa0(k)+LX(3)+aa0(k)*LX(4)+bb0(k)*LX(5);
        aa1(k+1)=aa1(k)+aa1(k)*LX(4)+bb1(k)*LX(5);
        aa2(k+1)=aa2(k)+aa2(k)*LX(4)+bb2(k)*LX(5);
        bb0(k+1)=bb0(k)+LX(6)+aa0(k)*LX(7)+bb0(k)*LX(8);
        bb1(k+1)=bb1(k)+aa1(k)*LX(7)+bb1(k)*LX(8);
        bb2(k+1)=bb2(k)+aa2(k)*LX(7)+bb2(k)*LX(8);
    end
    if max(Lo)<0.95
        continue;
    end
    
    XT=0;
    YT=0;
    GX=0;
    GY=0;
    for i=2:2*r+2
        for j=2:2*r+2
            gy=(Warp1(i,j+1)-Warp1(i,j-1))/2;
            gx=(Warp1(i+1,j)-Warp1(i-1,j))/2;
            XT=XT+wx(i)*gx^2;
            YT=YT+wy(j)*gy^2;
            GX=GX+gx^2;
            GY=GY+gy^2;
        end
    end
    
    XXT=XT/GX;
    YYT=YT/GY;
    WpointR(num_Bpoint,1)=XXT;
    WpointR(num_Bpoint,2)=YYT;
    
    
    BXT=aa0(k)+aa1(k)*XXT+aa2(k)*YYT;
    BYT=bb0(k)+bb1(k)*XXT+bb2(k)*YYT;
    BpointR(num_Bpoint,1)=BXT;
    BpointR(num_Bpoint,2)=BYT;
    
end

 %ind=find(BW);
 index=find(WpointR(:,1)==0);
  WpointR(index,:)=[];
  BpointR(index,:)=[];
  
% Bblkpos(:,1)=BpointR(:,2);
% Bblkpos(:,2)=BpointR(:,1);
% Wblkpos(:,1)=WpointR(:,2);
% Wblkpos(:,2)=WpointR(:,1);  
%   
% BW1=[Bblkpos Wblkpos];
% dlmwrite('E:\02cBeijing\72457\myfilename0824lms.txt', BW1, '\t');




