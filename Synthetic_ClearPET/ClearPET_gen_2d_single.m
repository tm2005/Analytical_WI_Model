function [rSector,module,crystal1,crystal2,layer,anglePhi] = ClearPET_gen_2d_single(angler,type)
% this gives us active sectors, modules, crystals and layers for one of two types of configuration

% type = 1 - 8*4 crystal on each side
% type = 2 - 8*2 crystal on each side

if type==1
    rSector = zeros(2,64);
    crystal1 = zeros(2,64);
    crystal2 = zeros(2,64);    
    layer = zeros(2,64); 
    
%     anglePhi = angler*ones(2,4096);
    k=1;
    for rs1=0:3
        for c21=0:7
            for l1=0:1
                rSector(1,k)=rs1;
                crystal2(1,k)=c21;
                layer(1,k)=l1;
                k=k+1;
                        
            end       
        end      
    end

    k=1;

    for rs2=10:13
        for c22=0:7   
            for l2=0:1
                rSector(2,k)=rs2;
                crystal2(2,k)=c22;
                layer(2,k)=l2;
                k=k+1;
            end         
        end   
    end
    

    ind = mod(rSector,2)==0 & mod(rSector,2)==0;
    crystal1(ind)=0;
    crystal1(ind)=0;
    ind = mod(rSector,2)==0 & mod(rSector,2)==1;
    crystal1(ind)=0;
    crystal1(ind)=4;
    ind = mod(rSector,2)==1 & mod(rSector,2)==0;
    crystal1(ind)=4;
    crystal1(ind)=0;
    ind = mod(rSector,2)==1 & mod(rSector,2)==1;
    crystal1(ind)=4;
    crystal1(ind)=4;

end

if type==2
    rSector = zeros(2,32);
    crystal1 = zeros(2,32);
    crystal2 = zeros(2,32);    
    layer = zeros(2,32); 

    k=1;

    for rs1=0:2:3
        for c21=0:7
            for l1=0:1
                rSector(1,k)=rs1;
                crystal2(1,k)=c21;
                layer(1,k)=l1; 
                k=k+1;         
            end    
        end 
    end

    k=1;
    
    for rs2=10:2:13     
        for c22=0:7
            for l2=0:1
                rSector(2,k)=rs2;
                crystal2(2,k)=c22;
                layer(2,k)=l2;
                k=k+1;
            end      
        end     
    end
    
    ind = mod(rSector,2)==0 & mod(rSector,2)==0;
    crystal1(ind)=0;
    crystal1(ind)=0;
    ind = mod(rSector,2)==0 & mod(rSector,2)==1;
    crystal1(ind)=0;
    crystal1(ind)=4;
    ind = mod(rSector,2)==1 & mod(rSector,2)==0;
    crystal1(ind)=4;
    crystal1(ind)=0;
    ind = mod(rSector,2)==1 & mod(rSector,2)==1;
    crystal1(ind)=4;
    crystal1(ind)=4;
    
end
anglePhi = angler*ones(size(layer));
module   = ones(size(layer));
end