close all
clear all
clc
%% Define PDF

R0 = 50; % 2R0 is the distance between two crystals
L0 = 10; % 2L0 is length of a crystal
h = 20; % distance from centar of rotation (origan) and line that connects centers of crystals

RR=R0*1.1;
LL=RR;

dr = 1e-1;
dz = 1e-1;
r = -RR:dr:RR;
z = -LL:dz:LL;

[R,Z] = meshgrid(r,z);

%% Expression is:

ind1 =  abs(R)<=RR & abs(Z-h)<=L0 & Z<=h & Z>=h-L0; 
ind2 =  abs(R)<=RR & abs(Z-h)<=L0 & Z>h & Z<h+L0; 
P = zeros(size(ind2));

P(ind1) = 1-abs(Z(ind1)-h)/L0;
P(ind2) = 1-abs(Z(ind2)-h)/L0;

%% Plot

f = figure;
f.Position = [100 100 800 800];
set(gcf, 'PaperSize', [9 9]);
surfl(R,Z,P),

colormap bone, shading interp, xlim([-RR RR]), ylim([-LL LL]), xlabel('r'), ylabel('z'), zlabel('P'), title('Non-normalized aproximation of the response of two crystals') ;
% saveas(gcf,'PDF2D_aprox.pdf')
%% Checking

Area = sum(P(:))*dr*dz; %must be close to 1
