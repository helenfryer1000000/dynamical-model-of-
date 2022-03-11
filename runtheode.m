function [TTshort,ZZshort,starty]=runtheode(tvecXi,startx,theparams)
% function [TTshort,ZZshort,TTy,ZZy,starty]=runtheode(tvecXi,startx,theparams)

global z  a k deltaL phi p B Bt beta

     z=theparams(1);
     a=theparams(2);     
     k=theparams(3);  
deltaL=theparams(4);  
   phi=theparams(5);
     p=theparams(6);
     B=theparams(7);
    Bt=theparams(8);
  beta=theparams(9);
     
[TTy,ZZy]=ode45(@dx,tvecXi,startx);

if(length(tvecXi)==2)
    TTshort=TTy([1,length(TTy)]);
    ZZshort=ZZy([1,length(TTy)],:);
else
    TTshort=TTy;
    ZZshort=ZZy;
end

sizeZZy=size(ZZy);
starty=ZZy(sizeZZy(1),:);