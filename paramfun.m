function ZZy = paramfun(paramsforguess,timething) %takes in the 

global  BSN BSP BtSN BtSP betaSN betaSP p n delta deltat Vindex pvec pvecSP tvecX tvecXSP fifLconv

justforformat=timething; %DON'T delete - required for formatting

zbig=paramsforguess(1);
a1=paramsforguess(2);
a2=paramsforguess(3);
k1=paramsforguess(4);
phi=paramsforguess(5);
k2=k1;

% zbig=paramsforguess(1);
% a1=paramsforguess(2);
% k1=paramsforguess(3);
% k2=paramsforguess(4);
% phi=paramsforguess(5);
% a2=a1;

dL2=p-a2+0.00053; dL1=dL2;
      zvec=[0    zbig   zbig  zbig   zbig zbig   zbig   zbig  0.3   0.3   zbig];
     avec=[a1    a1     a1    a1     a1   a1     a1     a1    a1    a1    a2  ];
     kvec=[k1    k1     k1    k1     k1   k1     k1     k1    k1    k1    k2 ];
deltaLvec=[dL1   dL1    dL1   dL1    dL1  dL1    dL1    dL1   dL1   dL1   dL2 ];
   phivec=[0     0      0     0      0    0      0      0     phi   phi   phi ];
   
  
  Istart=zeros(n+1,1); Istart(1)=1;
 Itstart=zeros(n+1,1);
  Vstart=zeros(n+1,1);
  Lstart=zeros(n+1,1);
  Tstart=BSN/delta;
 Ttstart=BtSN/deltat;
 
 TstartSP=BSP/delta;
 TtstartSP=BtSP/deltat;
 
LL=length(avec);
   TTcell=cell(LL,1);
   ZZcell=cell(LL,1);
startcell=cell(LL+1,1);
startcell{1}=[Istart;Itstart;Vstart;Lstart;Tstart;Ttstart];


TT=[];
ZZ=[];
for i=1:LL
    para=[zvec(i),avec(i),kvec(i),deltaLvec(i),phivec(i),pvec(i),BSN,BtSN,betaSN];
    [TTcell{i},ZZcell{i},startcell{i+1}]=runtheode(tvecX{i},startcell{i},para);
    TT=[TT;TTcell{i}(2:length(TTcell{i}))];
    ZZ=[ZZ;ZZcell{i}(2:length(TTcell{i}),:)];
end

     zvecSP=[0      zbig];
     avecSP=[a1     a2  ];     
     kvecSP=[k1     k2  ];
deltaLvecSP=[dL1    dL2 ];
   phivecSP=[phi    phi];
   
   TTcellSP=cell(2,1);
   ZZcellSP=cell(2,1);
startcellSP=cell(3,1);
startcellSP{1}=[Istart;Itstart;Vstart;Lstart;Tstart;Ttstart];  
   
TTSP=[];
ZZSP=[];
for i=1:2
    paraSP=[zvecSP(i),avecSP(i),kvecSP(i),deltaLvecSP(i),phivecSP(i),pvecSP(i),BSP,BtSP,betaSP];
    [TTcellSP{i},ZZcellSP{i},startcellSP{i+1}]=runtheode(tvecXSP{i},startcellSP{i},paraSP);
    TTSP=[TTSP;TTcellSP{i}(2:length(TTcellSP{i}))];
    ZZSP=[ZZSP;ZZcellSP{i}(2:length(TTcellSP{i}),:)];
end

lengthTT=length(TT);
lengthTTSP=length(TTSP);

 Vout=ZZ(:,Vindex);
  sumV=sum(Vout, 2);
Vdiv=divfun(Vout,lengthTT);

 VoutSP=ZZSP(:,Vindex);
  sumVSP=sum(VoutSP, 2);
 VdivSP=divfun(VoutSP,lengthTTSP);
  VdivSPdiff=VdivSP(1)-VdivSP(3);
 
ZZy=[log10(sumV/fifLconv);200*log10(sumVSP(3)/fifLconv);1000*Vdiv([6;13;26]);8000*VdivSPdiff];




