function dx=dx(t,x)

global delta deltat deltaI deltaIt deltaL pi z k p a c gamma
global Iindex Itindex Vindex Lindex Tindex Ttindex phi B Bt beta

 I=x(Iindex);
 It=x(Itindex);
 V=x(Vindex);
 L=x(Lindex);
 T=x(Tindex);
 Tt=x(Ttindex);

theta=thetafun(V);

 dI=(1-z)*(1-phi)*T*beta.*theta + gamma*It +a*L-deltaI*I;
dIt=(1-z)*(1-phi)*(1-k)*Tt*beta.*theta-gamma*It-deltaIt*It;
dV=pi*I-c*V;
dL=(1-z)*(1-phi)*k*Tt*beta.*theta+(p-a-deltaL)*L;
dT=B-(1-z)*(1-phi)*T*beta.*sum(V)-delta*T;
dTt=Bt-(1-z)*(1-phi)*Tt*beta.*sum(V)-deltat*Tt;

dx=[dI;dIt;dV;dL;dT;dTt];

end