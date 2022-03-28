close all
clear all
clc
%% Define PDF


R0 = 50;
L0 = 10;
Res =  sqrt(R0^2+L0^2)*1.1;
dr = 1e-1;
dz = 1e-1;

r = linspace(-Res,Res,1000);
z = linspace(-Res,Res,1000);

[R,Z] = meshgrid(r,z);

%% Expression is:

ind1 = abs(Z) <= L0./R0.*abs(R) & abs(R)<=R0+0*dr & abs(Z)<=L0;
ind2 = abs(Z) > L0./R0.*abs(R) & abs(Z)<=L0 & abs(R)<=R0+0*dr;

P = zeros(size(ind1));

P(ind1) = R0./(R0+abs(R(ind1))) ;
P(ind2) = R0^2./(R0^2-abs(R(ind2)).^2).*(L0-abs(Z(ind2)))./L0;

P=P/2/L0/R0;

%%

figure, imagesc(P), axis equal, colormap gray;
figure, surfl(P), colormap bone,shading interp;
%% Rotation

P_rot = zeros(size(P));
angle_step=0.1;
angle_range = 0:angle_step:360-angle_step;
for angle = angle_range
    P_rot = P_rot + imrotate(P,angle,'bilinear','crop');
end
P_rot = P_rot/length(angle_range);
figure, imagesc(P_rot), axis equal, colormap gray;
figure, surfl(P_rot), colormap bone,shading interp;

Area_rot = sum(P_rot(:))*dr*dz; %must be close to 1
%% finding P(r)

rr = sqrt(R.^2+Z.^2);
rr=rr(:);
[rr, ind] = sort(rr);
P_rot=P_rot(ind);
figure, plot(rr(:),P_rot(:));
%%
save('highres_numerical','rr','P_rot','R0','L0')