%x����Ҫ�˲���ͼ��,n��ģ���С(��n��n)  
function d = avg_filter(x,n)     
a(1:n,1:n)=1;   %a��n��nģ��,Ԫ��ȫ��1  
[height, width]=size(x);   %����ͼ����hightxwidth��,��hight>n,width>n  
x1=double(x);  
x2=x1;  
for i=1:height-n+1  
    for j=1:width-n+1  
        c=x1(i:i+(n-1),j:j+(n-1)).*a; %ȡ��x1�д�(i,j)��ʼ��n��n��Ԫ����ģ�����  
        s=sum(sum(c));                 %��c�����и�Ԫ��֮��  
        x2(i+(n-1)/2,j+(n-1)/2)=s/(n*n); %����ģ�������ĸ�Ԫ�صľ�ֵ����ģ������λ�õ�Ԫ��  
    end  
end
d= x2; 