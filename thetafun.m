function thetafun=thetafun(VV)

global n mu

VVstart=VV;
VVstart(length(VV))=[];
VVend=VV;
VVend(1)=[];

vvstart=0:(n-1);
vvstart=vvstart';
vvend=1:n;
vvend=vvend';

theta0=VV(1)*(1-mu*n);
thetarest=VVstart*mu.*(n-vvstart)+VVend.*(1-mu*(n-vvend));
thetafun=[theta0;thetarest];

