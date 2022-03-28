close all
clear all
clc
%%
const_list

load('nema30s.mat'); 
load('allzn');

X3D = zeros(256,256,48);
maxiter = 40;

%% From sectors, modules, crystals and layero to Cart. coordinates. (take a look at "CompToCartesian.m")

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
zn=zn -max(zn(:))/2; %centering

clearvars x y z
%% Slice by slice
for i = 1:48
%%  Take events only in allzn(i) plane
    ind1 = ( abs(zn(1,:) - allzn(i)) <1e-2);
    ind2 = ( abs(zn(2,:) - allzn(i)) <1e-2);
    ind = ind1 & ind2; % all indices on axial distance
    
    xi = xn(:,ind);
    yi = yn(:,ind);
    
    x=xi; %change of variables
    y=yi; %change of variables (this should be erased, but...)
%%  Preparing the WI 
    xm1 = x(1,:)';
    xm2 = x(2,:)';
    
    ym1 = y(1,:)';
    ym2 = y(2,:)';
    if ismember(i,[5;6;7;8;17;18;19;20;29;30;31;32;41;42;43;44])
        Ic =load("comp_256_type1_075.mat",'I1');
    else
        Ic =load("comp_256_type2_075.mat",'I1');
    end
    
    WI= Ic.I1;

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
    degstep=1;
    phideg =  -90+degstep:degstep:90;
    
    
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
    
%         figure(5), imagesc(X), colormap gray,axis equal;    % Current reconstruction
    
    end    
    X=X/max(X(:));
    X3D(:,:,i)=X;
    name="set1/img" + i + ".png";
    imwrite(X,name);
    i
end