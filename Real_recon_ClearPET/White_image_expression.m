close all
clear all
clc
%%
angler=0;
type =1;
angle0 = 0;
%% Geom
R0  = 67.8*3/4;
Rs1 = 72.8; %+5
Rs2 = 82.8; %+5+10
L0 = 1;

[rSector,module,crystal1,crystal2,layer,anglePhi] = ClearPET_gen_2d_single(0,type);
[xn,yn,~,~,~,~] = CompToCartesian(rSector, module, crystal1,crystal2,layer,anglePhi);

k=0;
for i = 1:size(xn,2)
    for j = 1:size(yn,2)
    k=k+1;
    R1(k)=sqrt(xn(1,j)^2+yn(1,j)^2);
    R2(k)=sqrt(xn(2,i)^2+yn(2,i)^2);
    Rij(k)=(sqrt( (xn(2,i)-xn(1,j)).^2 +  (yn(2,i)-yn(1,j)).^2 ))/2;
    hij(k)=abs( (yn(2,i)-yn(1,j))*xn(1,j) - (xn(2,i)-xn(1,j))*yn(1,j) )/(sqrt( (xn(2,i)-xn(1,j)).^2 +  (yn(2,i)-yn(1,j)).^2 ));
    end
end

Lij1 = L0*sqrt(1-hij.^2./R1.^2);
Lij2 = L0*sqrt(1-hij.^2./R2.^2);

Lij=mean([Lij1; Lij2]);

[P,r] = rot_trian2(hij,Rij,Lij,50000);

figure, plot(r,P);
%%

clear S
Nout = 256;
k0 = 63;

v = linspace(0,R0-R0/128*sqrt(2),(Nout+2)/2);
d0 = mean(diff(v));
w = [-v(end-1:-1:2),v];

dv = mean(diff(v));

wsx = w;
wsy = w;
vs = linspace(0,dv/2*(k0-1)/(k0+1),(k0+1)/2);
ws = [-vs(end:-1:2),vs];

[Xs,Ys] = meshgrid(ws, ws);

Iout= zeros(length(w),length(w));
for i=1:1:length(w)
    for j=1:1:length(w)
        Xl = wsx(i)*ones(k0,k0);
        Yl = wsy(j)*ones(k0,k0);
        Xl = Xl+Xs;
        Yl = Yl+Ys;
        Rl = sqrt(Yl.^2+Xl.^2);
        Il = interp1(r,P,Rl(:),'linear');
        Iout(i,j)=mean(Il(:));

    end
end


I1=Iout/max(Iout(:));
figure, imagesc(I1), colormap gray, axis equal;
% save('comp_256_type1_075','I1');
