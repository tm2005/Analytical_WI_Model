close all
clear all
clc
%% 
% This script simulates one intersection of ClearPET
% Rs1 and Rs2 are radii of to layer1 and layer2
% 2L0 is crystals length
% R0 is ROI
% type1 =  8 active sectors (8*2*4 ctystals on each side of a ring)
% type2 =  4 active sectors (8*2*2 ctystals on each side of a ring)

%% Geom
R0  = 67.8*3/4;
Rs1 = 72.8; %+5
Rs2 = 82.8; %+5+10
L0 = 1;
type = 2;  

angstep=1/10; %step of mechanical rotation

xm1=[];
xm2=[];
ym1=[];
ym2=[];


for angle0=0:angstep:360-angstep
    [rSector,module,crystal1,crystal2,layer,anglePhi] = ClearPET_gen_2d_single(angle0*10,type);
    [xn,yn,xun,yun,xln,yln] = CompToCartesian(rSector, module, crystal1,crystal2,layer,anglePhi);

    %% CHOOSE A PHANTOM:
    %% 1) Uniform image (for white (compensations) image)
    
%     Nevents = 50;
%     [x0,y0]=generate_radnom_events_for_circular_white_image((2.6*(type-1)+1.3)*Nevents,30);

    %% 2) A part of NEMA phantom with 5 small circles

%     N0=200*3.5;
%     Radius_large=33.5;
%     [x0,y0]=nema_nu4_five(N0,Radius_large/2/2,pi/12,[1;2;3;4;5]/2,0,0);

    %% 3) A part of NEMA phantom with 2 empty circles (holes)

    N0=260;
    Radius_large=30;
    Radius_hole=8;
    [x0,y0] = nema_nu4_holes(N0,Radius_large/2,pi/9,Radius_large/2/2,[Radius_hole;Radius_hole]/2,0,0);

    %%

    phiwi = rand(length(x0),1)*pi - pi/2; 
    kwi = tan(phiwi); % random angles for rays
    
    for  i = 1: length(x0)
        for j = 1:size(yun,2)

        D = -kwi(i)*(xun(1,j)-xln(1,j)) + yun(1,j)-yln(1,j);
        s(1) = ( (-kwi(i)*x0(i)+y0(i)-yln(1,j) )*(xun(1,j)-xln(1,j)) + ( yun(1,j)-yln(1,j) )*xln(1,j) )/D;
        s(2) = ( (kwi(i)*(xln(1,j)-x0(i))+y0(i))*(yun(1,j)-yln(1,j)) - kwi(i)*yln(1,j)*(xun(1,j)-xln(1,j)) )/D;

        for k = 1:size(yun,2)

            D = -kwi(i)*(xun(2,k)-xln(2,k)) + yun(2,k)-yln(2,k);
            p(1) = ( (-kwi(i)*x0(i)+y0(i)-yln(2,k) )*(xun(2,k)-xln(2,k)) + ( yun(2,k)-yln(2,k) )*xln(2,k) )/D;
            p(2) = ( (kwi(i)*(xln(2,k)-x0(i))+y0(i))*(yun(2,k)-yln(2,k)) - kwi(i)*yln(2,k)*(xun(2,k)-xln(2,k)) )/D;
            
            if s(1)<=max(xun(1,j),xln(1,j)) && min(xun(1,j),xln(1,j))<=s(1) && s(2)<=max(yun(1,j),yln(1,j)) && min(yun(1,j),yln(1,j))<=s(2)
                if p(1)<=max(xun(2,k),xln(2,k)) && min(xun(2,k),xln(2,k))<=p(1) && p(2)<=max(yun(2,k),yln(2,k)) && min(yun(2,k),yln(2,k))<=p(2)
    
         	    xm1 = [xm1; xn(1,j)];
                xm2 = [xm2; xn(2,k)];
                ym1 = [ym1; yn(1,j)];
                ym2 = [ym2; yn(2,k)];
                end
            end
        end
        end
    end
    angle0 %counter
end

save('measurementxyz');