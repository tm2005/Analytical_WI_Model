close all
clear all
clc
%%
angler=0;
type = 2;
% angle0 = 0;
%% Geom
R0  = 67.8*3/4;
Rs1 = 72.8; %+5
Rs2 = 82.8; %+5+10
L0 = 1;
I=zeros(128,128);

parfor angle0=0:1800-1
    [rSector,module,crystal1,crystal2,layer,anglePhi] = ClearPET_gen_2d_single(angle0,type);
    [xn,yn,xun,yun,xln,yln] = CompToCartesian(rSector, module, crystal1,crystal2,layer,anglePhi);
    % figure, scatter(xn(1,:),yn(1,:)), hold on, scatter(xn(2,:),yn(2,:)),  axis equal
    % figure, scatter(xun(1,:),yun(1,:)), hold on, scatter(xln(1,:),yln(1,:)), scatter(xn(1,:),yn(1,:)), axis equal
    
    %% 
    %Gen events
    Nevents = 18000;
    [x0,y0]=generate_radnom_events_for_circular_white_image(Nevents,R0);
    phiwi = rand(Nevents,1)*pi - pi/2;
    kwi = tan(phiwi);
    l=0;
    xm1=[];
    xm2=[];
    ym1=[];
    ym2=[];
    for  i = 1: Nevents
        for j = 1:size(yun,2)
        s = zeros(2,1);

        D = -kwi(i)*(xun(1,j)-xln(1,j)) + yun(1,j)-yln(1,j);

        s(1) = ( (-kwi(i)*x0(i)+y0(i)-yln(1,j) )*(xun(1,j)-xln(1,j)) + ( yun(1,j)-yln(1,j) )*xln(1,j) )/D;
        s(2) = ( (kwi(i)*(xln(1,j)-x0(i))+y0(i))*(yun(1,j)-yln(1,j)) - kwi(i)*yln(1,j)*(xun(1,j)-xln(1,j)) )/D;
        for k = 1:size(yun,2)
            p = zeros(2,1);

            D = -kwi(i)*(xun(2,k)-xln(2,k)) + yun(2,k)-yln(2,k);

            p(1) = ( (-kwi(i)*x0(i)+y0(i)-yln(2,k) )*(xun(2,k)-xln(2,k)) + ( yun(2,k)-yln(2,k) )*xln(2,k) )/D;
            p(2) = ( (kwi(i)*(xln(2,k)-x0(i))+y0(i))*(yun(2,k)-yln(2,k)) - kwi(i)*yln(2,k)*(xun(2,k)-xln(2,k)) )/D;
            
            if s(1)<=max(xun(1,j),xln(1,j)) && min(xun(1,j),xln(1,j))<=s(1) && s(2)<=max(yun(1,j),yln(1,j)) && min(yun(1,j),yln(1,j))<=s(2)
                if p(1)<=max(xun(2,k),xln(2,k)) && min(xun(2,k),xln(2,k))<=p(1) && p(2)<=max(yun(2,k),yln(2,k)) && min(yun(2,k),yln(2,k))<=p(2)
                l=l+1;

                a = rand;
                T1 = a*[xun(1,j);yun(1,j)]+(1-a)*[xln(1,j);yln(1,j)];
                a = rand;
                T2 = a*[xun(2,k);yun(2,k)]+(1-a)*[xln(2,k);yln(2,k)];
    
                xm1 = [xm1;T1(1)];
                xm2 = [xm2;T2(1)];
                ym1 = [ym1;T1(2)];
                ym2 = [ym2;T2(2)];

                end
            end
        end
        end
    end
    %% Backprojection

    xi1=xm1;
    xi2=xm2;
    yi1=ym1;
    yi2=ym2;
    
    K = (yi2-yi1) ./ (xi2-xi1);
    L = yi1 - K .* xi1;
    
    % r and phi for sinogram from line
    r2d = L ./ sqrt(1 + K.^2);
    phi = atan((yi2-yi1) ./ (xi2-xi1));
    
    Sc = [r2d phi];
  
    Sg = hist3(Sc,{-R0*sqrt(2):R0/91*sqrt(2):R0*sqrt(2), -pi/2:pi/360:pi/2});
    I0 = iradon(Sg, 180*[-pi/2:pi/360:pi/2]/pi,'None'); %FBP
    I = I0 + I;
    display(angle0)
end

figure, imagesc(I), colormap gray, axis equal
%%
I1=I/max(I(:));
save('comp_MC_32_128_R075_xyz','I1');
