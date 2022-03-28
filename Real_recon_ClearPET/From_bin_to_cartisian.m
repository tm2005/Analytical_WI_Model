close all
clear all
clc
%%
load('allzn');
const_list % all geometrical constants

%%
% Files nema10s.mat, nema10s.mat and Norm140322.mat containt information about all events in form of sectors, modules, crystals, layers and angle.
% This script turns such information to Cartesian coordinates. 
% Choose a file:
% load('nema30s.mat'); 
load('Norm140322.mat'); 

% choose intersection (from 1 to 48) (README.txt) 
loc = 3;

%% 
x=zeros(size(rSector));
y=zeros(size(rSector));
z=zeros(size(rSector));

phi =  double(anglePhi)*azimuthalStep*pi/180 - pi/2;

x = x + (radius + crystalDepth/2);
y = y + (double(crystal2) - 7/2)*crystalPitch;
z = z + (double(crystal1))*crystalPitch;

ind = find(layer==1);
x(ind) = x(ind) + crystalDepth;

phi = phi + double(rSector)*sectorPitch*pi/180;

xn = x;

xn = x.*cos(phi)-y.*sin(phi);
yn = y.*cos(phi)+x.*sin(phi);

zn = z + double(module)*(modulePitch);
zn = zn + shift*mod(double(rSector+1),2);
zn=zn -max(zn(:))/2;
%% Choosing intersection

clearvars x y z
ind1 = ( abs(zn(1,:) - allzn(loc)) <1e-2);
ind2 = ( abs(zn(2,:) - allzn(loc)) <1e-2);
ind = ind1 & ind2; % all indices on axial distance

xi = xn(:,ind);
yi = yn(:,ind);

x=xi; %change of variables
y=yi; %change of variables (this should be erased, but...)

save('measure_WI_3','x','y')
