%Task 1.3
%syms eta
Pri=0.86;
C=0.2;%
%Prt=@(eta) 
a=1;
b=1e4;
n=1e7;
eta=linspace(a,b,n);
dx=(b-a)/n;
I(1)=0;
for i=1:length(eta)
    if i==1
        j=2;
    elseif i>1
        j=i;
    end
    
I(i)=I(j-1)+dx*1/3*eta(i)^(-2/3)*(1+eta(i)*((.5*Pri.^-1+eta(i)*C*(Pri.^(-.5))-C^2*eta(i)^2*(1-exp(-1/C*Pri.^(.5)*eta(i)^(-1))))^-1)^-1)^-1;
end
%I=quadgk(Prt,1e-3,inf);
%%
plot(eta,I)
%%
%Task 1.3
%syms eta
% Pri=0.86;
% C=0.2;%
% Prt=@(eta) (.5*Pri.^-1+eta*C*(Pri^(-.5))-C^2*eta.^2*(1-exp(-1/C*Pri^(.5)*eta.^(-1))));
I=quadgk(@Prt,1e-3,inf);
% a=1;
% b=1e4;
% n=1e7;
% eta=linspace(a,b,n);
% dx=(b-a)/n;
% I(1)=0;
% for i=1:length(eta)
%     if i==1
%         j=2;
%     elseif i>1
%         j=i;
%     end
%     
% I(i)=I(j-1)+dx*
% end
%I=quadgk(Prt,1e-3,inf);
%%
r=5.1e-4;
C16=r^(1/3)/I/(245/8)^(7/24);
%%
Pr=[10 100 1000];
Re=linspace(1e2,1e5);
for i=1:3
    hk1=.056*Pr(i).^(1/3)*Re.^.2;
    hk2=C16*Pr(i).^(1/3)*Re.^(7/24);
%     ssres(i)=sum((hk1-hk2).^2);
%     sstot(i)=sum(hk1-mean(hk1));
%     R2(i)=1-ssres(i)/sstot(i)
    figure;
    plot(Re,hk1,'--',Re,hk2);
    xlabel('Re');ylabel('LHS');%title('Pr=%d',Pr(i));
    legend('Colburn','Eq (1.13)');
end
%%
alpha=zeros(3,4); RJalPrl=alpha; nd=alpha; 
mginfmgi=alpha; mgi=alpha; PviPinf=alpha; 
Ti=alpha; Tw=alpha; vi=alpha;
F0=[.4 1.25 5.5];
mginf=linspace(.004,.08,4);
for i=1:3
    for j=1:4
[alpha(i,j) RJalPrl(i,j) nd(i,j) mginfmgi(i,j) mgi(i,j) PviPinf(i,j) Ti(i,j) Tw(i,j) vi(i,j) NuRe(i,j)]=fnp(F0(i),mginf(j));
    end
%     figure;
%     plot(mginf,Tw(i,:));
%         figure;
%     plot(mginf,Ti(i,:));
%         figure;
%     plot(mginf,mginfmgi(i,:));
%         figure;
%     plot(mginf,mgi(i,:));
%         figure;
%     plot(mginf,NuRe(i,:));
end
%%
figure;
plot(mginf,RJalPrl(1,:),mginf,RJalPrl(2,:),mginf,RJalPrl(3,:))
Ti=Ti-273;
plot(mginf,Ti(1,:),mginf,Ti(2,:),mginf,Ti(3,:))





























