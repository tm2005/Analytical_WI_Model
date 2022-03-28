close all
clear all
clc
%% This is MLEM-like reconstruction

%% Geometry
R0  = 67.8*3/4;
Rs1 = 72.8; %+5
Rs2 = 82.8; %+5+10
%%

maxiter=40; % total number of iterations
%%
% Load measurement
load('measure_NEMA_18.mat') %comp_type1
% load('measure_NEMA_36.mat') %comp_type2

%%
xm1 = x(1,:)';
xm2 = x(2,:)';

ym1 = y(1,:)';
ym2 = y(2,:)';

%%
% Load white iamge
Ic =load("comp_256_type1_075.mat",'I1');
% Ic =load("comp_256_type2_075.mat",'I1');
% 
WI= Ic.I1;

figure, imagesc(WI), colormap gray, axis equal %white image
n=size(WI,2);
WI=double(WI);
%%

v = linspace(0,(R0),(256/2+1));
w = [-v(end-1:-1:2),v];
[X,Y]=meshgrid(w,w);
R=sqrt(X.^2+Y.^2);
ind0 = R<R0;

WI=WI./max(WI(:));

%%
% Init guess:

X = 1./WI; % this one
% X = ones(size(WI)); %or this one

X(~ind0)=0;
%% Define image space

v = linspace(0,(R0)*(sqrt(2)),(367/2+1));
w = [-v(end:-1:2),v];
dw = mean(diff(w));
ind0 =  abs(WI)>1e-16;
we = [w-dw/2];
degstep=0.5;
phideg =  -90+degstep:degstep:90;

tic
for k=1:maxiter % Iterations
%%  Dithering
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

%%  Create sinogram Sm
    K = (y2-y1) ./ (x2-x1);
    L = y1 - K .* x1;

    r2 = L ./ sqrt(1 + K.^2);
    phi = atan((y2-y1)./ (x2-x1));
    
    Sc = [r2 phi]; 
    Sm = hist3(Sc,'Edges',{we, (phideg-1/2*degstep)*pi/180});

    Sm(:,1)=0;
    Sm(:,end)=0;
%%  MLEM-like recon

    den = radon(X,phideg);

    ind = Sm< 1e-16 & den < 1e-16;
    S = Sm./den;
    S(ind) = 0;

    Xs = iradon(S,phideg,'nearest','None');

   	Xs = wkeep(Xs,size(X));
    X(ind0) = X(ind0).*Xs(ind0)./WI(ind0);

    figure(5), imagesc(X), colormap gray,axis equal;    % Current reconstruction

end
toc
X=X/max(X(:));

% imwrite(X,'Ideal_recon_type1_five.png')

