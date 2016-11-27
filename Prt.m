function [Prtt] = Prt(eta)
Pri=.86;
C=.2;

Prtt=1/3*eta.^(-2/3).*(1+eta.*((.5*Pri.^-1+eta.*C*(Pri.^(-.5))-(C^2*eta.^2).*(1-exp(-1/C*Pri.^(.5)*eta.^(-1)))).^-1).^-1).^-1;


end

