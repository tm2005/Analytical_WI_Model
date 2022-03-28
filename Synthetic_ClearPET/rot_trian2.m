function [P,r] = rot_trian2(hij,Rij,Lij,Nr)
% Summing over all crystal-to-crystal responses (triangular rot_trian_1d_Rij)
    rs = max(Rij)/(Nr-1);
    r = 0:rs:2*max(Rij);
    P=zeros(size(r));
    L = length(hij);
    for i = 1:L
        h = hij(i);
        R = Rij(i);
        L = Lij(i);
        P = P + rot_trian_1d_Rij(r,h,R,L);
    end
    P = P/length(hij);
end