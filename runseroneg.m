clear all
clc

%lsqnonlin more general fitting?

global R0 delta deltat deltaI deltaIt pi c tvec tvecSP pvec pvecSP n tvecX tvecXSP p fifLconv gamma
global Iindex Itindex Vindex Lindex Tindex Ttindex betaSN betaSP BSN BSP BtSN BtSP mu
n=100;

  Iindex=1:(n+1);
 Itindex=(1:(n+1))+ Iindex(length(Iindex));
  Vindex=(1:(n+1))+Itindex(length(Itindex));
  Lindex=(1:(n+1))+ Vindex(length(Vindex));
 Tindex=Lindex(length(Lindex))+1;
Ttindex=Tindex(length(Tindex))+1;

mu=10^(-5);
fifLconv=7500;
R0=4.6;
delta=0.15;
deltaIt=0.01;
deltat=0.01;
c =23;
deltaI=0.7;
gamma=0.05; %need to check that this matches up with latest refs etc
pi=6000;
p=0.02;
BSN=310*fifLconv; %seronegative (activated)
BSP=700*fifLconv; %seropositive (activated)
BtSN=BSN/100; 
BtSP=BSP/100;
betaSN=R0*delta*deltaI*c/BSN/pi; % seronegative(resting)
betaSP=R0*delta*deltaI*c/BSP/pi; % seropositive(resting)

tvecX{1} =[-120 0];
tvecX{2} =[0    15   40   75   109  146];
tvecX{3} =[146  214];
tvecX{4} =[214  256];
tvecX{5} =[256  279  306  350  469];
tvecX{6} =[469  560];
tvecX{7} =[560  637];
tvecX{8} =[637  706  908  1022 1154 1202 1218 1274 1355 1468 1621 1716];
tvecX{9} =[1716 1883];
tvecX{10}=[1883 1926];
tvecX{11}=[1926 1951];

tvecXSP{1} =[-1500 -183  0];
tvecXSP{2} =[0     2520];
  
tvec=[-120  0     146   214    256  469    560    637  1716   1883  1926 3000];
pvec=[p     p     p     p      p    p      p      p    p      p     p  ];
tvecSP=[-1500  0      2520];   
pvecSP=[p     p]; 


%% this changes activation
zsmallG=0.6055;
zbigG=0.7856;
a1G=0.004598; 
a2G=1.972e-05;
k1G=0.0016; 
phiG=0.63;

paramguess=[zbigG    a1G    a2G       k1G     phiG];
  lowparam=[0.78     0.004  0.000015  0.0015  0.625]; 
   upparam=[0.79     0.005  0.000025  0.0017  0.635];
%    

%%%this changs k (production of latent cells)
% zsmallG=0.6055;
% zbigG=0.7856;
% a1G=1.972e-05;
% k1G=0.0016;
% k2G=0.000016;
% phiG=0.63;
% paramguess=[zbigG    a1G     k1G       k2G       phiG];
%   lowparam=[0.78     0.003   0.014     0.000014   0.63]; 
%    upparam=[0.79     0.006   0.018     0.000018   0.64];

xxxxx=[0,5.5,15.5,16.5,19.5,2,-1,7.5,17,13,9.5,8,14,1,8.5,8,4,19,2.5,20,18,16,0,18,18.5,7.4,7.5,2.5]/22.5;
yyyyy=[5,5,  4,   4,   4,   4,4, 5,  3, 3, 3,  3,4, 4,3,  3,3,2, 3,  2, 2, 2, 3,2, 2,   4,  4,  2];
VLdataclosed=xxxxx+yyyyy; %from the siemieniuk
VLdaysclosed=[0    15   40   75   109  146  214  256  279  306  350  469  560  637  706  908  1022 1154 1202 1218 1274 1355 1468 1621 1716 1883 1926 1951]; % from joe grove
VLdataopen=[2.1111 1.60 1.6  1.60 1.6  1.60 1.6  1.60 1.6  1.60 1.6  1.60];
VLdaysopen=[1951 2035 2140 2276 2349 2412 2469 2608 2678 2731 2819 2982];

palmerdayssmall=7*[0  0.5       1         1.5       2         2.5       3         3.5       4     ];
  palmerdatasmall=[5  4.4516    3.8065    3.4194    3.2258    3.1290    3.0968    3.0645    3.0000];

palmerdaysbig=7*[60     90     120    150    180    210    240    270    300    330    360];
palmerdatabig = [0.7500 0.5962 0.4615 0.3654 0.2885 0.2692 0.2308 0.2115 0.1923 0.1923 0.1923];
    
  nehertime=[0:500:3000];
neherSYNdiv=[0.0007    0.0024    0.0036    0.0048    0.0062    0.0073    0.0086];

meands=[0.0019    0.0015    0.0032];

DATAPOS=0.1923;
TIMEPOS=360*7;

THEDATA=[VLdataclosed';200*DATAPOS;1000*meands';0];
THETIME=[VLdaysclosed';TIMEPOS;[146;560;1883];0];

% THETIME=[VLdataclosed';VLdataopen';meands']; % make sure these days are correct!
% THEDATA=[VLdaysclosed';VLdaysopen';[146;560;1883]]; 

dat1index=[1:length(VLdataclosed)];
dat2index=length(VLdataclosed)+1;
dat3index=length(VLdataclosed)+1+[1:3];
dat4index=length(THEDATA);

[pbest,presnorm,presidual,exitflag,output] = lsqcurvefitH(@paramfun,paramguess,THETIME,THEDATA,lowparam,upparam);
odes = presidual + THEDATA;
fprintf(     'New parameters: %f, %f, %f, %f, %f, %f',pbest(1:6));
fprintf('Original parameters: %f, %f, %f, %f, %f, %f',paramguess);

save('newoutput')

figure(1)
hold on; plot(THETIME(dat1index),THEDATA(dat1index),'ro-','markerfacecolor','r') % data
hold on; plot(THETIME(dat1index),odes(dat1index),'rx:') % best fit
xlabel('Time since onset of therapy (days)')
ylabel('log_1_0(plasma VL/ml)')

figure(2)
hold on; plot(THETIME(dat2index),THEDATA(dat2index),'bs')
hold on; plot(THETIME(dat2index),odes(dat2index),'bx')
xlabel('Time since onset of therapy (days)')
ylabel('log_1_0(plasma VL/ml)')

figure(3)
hold on; plot(THETIME(dat3index),THEDATA(dat3index),'bo-')
hold on; plot(THETIME(dat3index),odes(dat3index),'bx:')
legend('data','best fit')
xlabel('time(days)')
ylabel('divergence')

D1=THEDATA(dat1index)-odes(dat1index)
D2=THEDATA(dat2index)-odes(dat2index)
D3=THEDATA(dat3index)-odes(dat3index)
D4=THEDATA(dat4index)-odes(dat4index)








