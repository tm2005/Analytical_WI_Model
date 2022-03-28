close all
clear all
clc
%%
% Load measurement
% load('Measurement_type1_twoholes_075_30min.mat')
% load('Measurement_type1_fivespikes_075_30min.mat')
load('Measurement_type2_twoholes_075_30min.mat')
% load('Measurement_type2_fivespikes_075_30min.mat')

%% Dithering
R1 = sqrt(xm1.^2+ym1.^2);

Rs = R1;
Rs(Rs<Rs2) = Rs1;
Rs(Rs>=Rs2) = Rs2;

phi1 = atan2(ym1,xm1);

D1 = sqrt(R1.^2+1-2*R1.*cos(asin(Rs./R1)));
dang1 = asin( sin(asin(Rs./R1))./D1);
D2 = sqrt(R1.^2+1-2*R1.*cos(pi-asin(Rs./R1)));
dang2 = asin( sin(pi-asin(Rs./R1))./D2)   ;
angle0 = (dang1+dang2)/2;

phid1 = angle0.*(2*rand(size(R1)) -1 );
phin1 = phi1+phid1;
R1n = Rs./cos(acos(Rs./R1)+phid1);

x1 = R1n.*cos(phin1);
y1 = R1n.*sin(phin1);

R2 = sqrt(xm2.^2+ym2.^2);

Rs = R2;
Rs(Rs<Rs2) = Rs1;
Rs(Rs>=Rs2) = Rs2;

phi2 = atan2(ym2,xm2);

D1 = sqrt(R2.^2+1-2*R2.*cos(asin(Rs./R2)));
dang1 = asin( sin(asin(Rs./R2))./D1);
D2 = sqrt(R2.^2+1-2*R2.*cos(pi-asin(Rs./R2)));
dang2 = asin( sin(pi-asin(Rs./R2))./D2)   ;
angle0 = (dang1+dang2)/2;

phid1 = angle0.*(2*rand(size(R2)) -1 );
phin2 = phi2+phid1;
R2n = Rs./cos(acos(Rs./R2)+phid1);

x2 = R2n.*cos(phin2);
y2 = R2n.*sin(phin2);

%% Sinogram
K = (y2-y1) ./ (x2-x1);
L = y1 - K .* x1;

r2 = L ./ sqrt(1 + K.^2);
phi = atan((y2-y1)./ (x2-x1));

v = linspace(0,(R0)*(sqrt(2)),(367/2+1));
w = [-v(end:-1:2),v];
dw = mean(diff(w));
we = w-dw/2;
degstep=0.1;
phideg =  -90+degstep:degstep:90;

Sc = [r2 phi]; 
Sm = hist3(Sc,'Edges',{we, (phideg-1/2*degstep)*pi/180});

figure, imagesc(Sm), colormap gray;

%% Reconstruction

Ibp = iradon(Sm,(phideg-1/2*degstep),'Shepp-Logan');
Ibp =Ibp/max(Ibp(:));

figure, imagesc(Ibp), colormap gray, axis equal;