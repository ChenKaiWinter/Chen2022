clear all;
close all;
clc
N=8;%阵元个数
f=300000000;%信号频率
c=300000000;
lambda=c/f;%波长为10
d=0.5*lambda;%单元间隔是1/2波长,就是5
theta=linspace(-90,90,181); %从-90到90，等间隔181个点
theta0=0; %预先设定的来波方向
k = 1;
p = zeros(3,N);%p是一个3*8的double类型的零矩阵
for iii = 0 : N-1
          p(:,k) = [0,iii*d,0]';
          k = k+1;
    
end   %阵元的位置
%p=[0  0    0   0    0    0    0    0
%   0  5    10  15   20   25   30   35
%   0  0    0   0    0    0    0    0]
figure(1);%画出线阵的空间分布图
plot3(p(1,:),p(2,:),p(3,:),'ko');  %  k:black;  o:circle;
hold on;
xlabel('\it x');
ylabel('\it y');
zlabel('\it z');

k0=-2*pi/lambda*[cos(theta0*pi/180),sin(theta0*pi/180),0].';% [cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)],phi=2*pi
%k0=[-pi/5
%     0
%     0]
w0=exp(1i*k0.'*p);  % 1i:复数
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
xlabel('方位角(deg)'), ylabel('阵列方向图(dB)')
axis([-90 90 -60 0]);%可为x轴和y轴设置一个极限范围