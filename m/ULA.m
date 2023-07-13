clear all;
close all;
clc
N=8;%��Ԫ����
f=300000000;%�ź�Ƶ��
c=300000000;
lambda=c/f;%����Ϊ10
d=0.5*lambda;%��Ԫ�����1/2����,����5
theta=linspace(-90,90,181); %��-90��90���ȼ��181����
theta0=0; %Ԥ���趨����������
k = 1;
p = zeros(3,N);%p��һ��3*8��double���͵������
for iii = 0 : N-1
          p(:,k) = [0,iii*d,0]';
          k = k+1;
    
end   %��Ԫ��λ��
%p=[0  0    0   0    0    0    0    0
%   0  5    10  15   20   25   30   35
%   0  0    0   0    0    0    0    0]
figure(1);%��������Ŀռ�ֲ�ͼ
plot3(p(1,:),p(2,:),p(3,:),'ko');  %  k:black;  o:circle;
hold on;
xlabel('\it x');
ylabel('\it y');
zlabel('\it z');

k0=-2*pi/lambda*[cos(theta0*pi/180),sin(theta0*pi/180),0].';% [cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)],phi=2*pi
%k0=[-pi/5
%     0
%     0]
w0=exp(1i*k0.'*p);  % 1i:����
for ii=1:length(theta)
   
        k=-2*pi/lambda*[cos(theta(ii)*pi/180),sin(theta(ii)*pi/180),0].';
        v=exp(1i*k.'*p);
        b(ii,:)=w0*v';
end
B=20*log10(abs(b)/max(max(abs(b))));

%c=abs(b);
%e=min(c);

figure(2);
plot(theta,B,'k-')
grid on,hold on
xlabel('��λ��(deg)'), ylabel('���з���ͼ(dB)')
axis([-90 90 -60 0]);%��Ϊx���y������һ�����޷�Χ