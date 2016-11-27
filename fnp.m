function [alpha RJalPrl nd mginfmgi mgi PviPinf Ti Tw vi NuRe]=fnp(F0,mginf)
global a b c d;
%F0=3;mginf=.02;
C21=.5;C22=.5;C23=2;C24=1;C25=2;C26=2;
pl=958;mul=2.84e-4;uinf=1.2;
kl=.682;cpl=4211;Prl=1.74;
pm=.598;num=1.21e-5;
Sc=.55;R=201.7;
Mg=28.97;Mv=18;
Pinf=101e3;Tinf=373;
alpha=.5;
E=1;
while E>1e-5
[T,Y]=ode45(@myfunc,[0 5],[F0 0 alpha 0 0]);
E=abs(Y(length(T),2)-1);
alpha102=1.02*alpha;
[T1,Y1]=ode45(@myfunc,[0 5],[F0 0 alpha102 0 0]);
E1=abs(Y1(length(T),2)-1);
alpha=alpha-E/(E1-E)*(alpha102-alpha);
end
xx=.5;
phip0=1/Y(length(T),5);%F50=y1(nx);
RJalPrl=((F0)^3/(C21*C24^2*alpha))^.5;
nd=(F0/(C21*alpha))^.5;
mginfmgi=1+(Sc*F0)/(C23*phip0);
mgi=mginf/mginfmgi;
PviPinf=((1-mgi)/(1-mgi*(1-Mv/Mg)));
Pvi=Pinf*PviPinf;
Ti=1./(log(Pvi)-b)*a;
hlv=Ti*c+d;
Jal=RJalPrl/R*Prl;
TiTw=Jal*hlv/cpl;
Tw=Ti-TiTw;
NuRe=((Ti-Tw)/(Tinf-Tw))/nd;
vi=-.5*(uinf*num/pm/xx)^.5*F0;
end

