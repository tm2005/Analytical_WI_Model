function [xn,yn,xun,yun,xln,yln] = CompToCartesian(rSector, module, crystal1,crystal2,layer,anglePhi)
% (xn,yn) center of a crystal
% (xun,yun) first edge of a crystal
% (xun,yun) second edge of a crystal

%% All constants in mm and deg

geometrical_design_type = 1;
ring_diameter = 135.6;
rsector_axial_pitch = 0;
rsector_azimuthal_pitch = 18 ;
rsector_tangential_size = 19.1;
rsector_axial_size = 111.1;
module_axial_size = 19.1;
module_tangential_size = 19.1; 
module_axial_pitch = 27.6;
module_tangential_pitch = 0;
submodule_axial_size = 19.1;
submodule_tangential_size = 19.1;
submodule_axial_pitch = 0;
submodule_tangential_pitch = 0;
crystal_axial_size = 2;
crystal_tangential_size = 2; 
crystal_radial_size = 20 ;
crystal_axial_pitch = 2.3;
crystal_tangential_pitch = 2.3;
layer0_radial_size = 10;
in_layer0_interaction_length = 5;
layer1_radial_size = 10;
in_layer1_interaction_length = 5;
clock_time_step = 390e-12;
azimuthal_step = -0.1;
z_shift_sector_0_mod_2 = 9.2;
axial_step = 0;

%% More constants

numSectors = 20;%	///< Number of sectors (=cassettes) in the scanner
numModules = 4;%	///< Number of modules (=PMTs) per sector 
numSubMods = 1;%	///< Number of submodules, which are not used in ClearPET systems
numPixelsAxial = 8;%	///< Number of pixels per module in axial direction
numPixelsTangential = 8;%	///< Number of pixels per module in tagential direction
numPixels = numPixelsAxial*numPixelsTangential;%	///< Total number of pixels per module
numLayers = 2;%	///< Number of layers per pixel
numRings = 48;%	///< total number of rings per scanner
numPixelsTotalTangential = numSectors*numPixelsTangential;


radius        = 0.5 * ring_diameter;
crystalSize   = crystal_tangential_size;
crystalDepth  = crystal_radial_size;
crystalDepth  = crystalDepth/2;
crystalPitch  = crystal_axial_pitch;
modulePitch   = module_axial_pitch;
sectorPitch   = rsector_azimuthal_pitch;
azimuthalStep = azimuthal_step;
shift         = z_shift_sector_0_mod_2;


%% Main part
x=zeros(size(rSector));
y=zeros(size(rSector));

phi =  double(anglePhi)*azimuthalStep*pi/180;

x = x + (radius + crystalDepth/2);

yo =  y + (double(crystal2) - 7/2)*crystalPitch;
yu =  y + (double(crystal2) - 7/2)*crystalPitch + 1;
yl =  y + (double(crystal2) - 7/2)*crystalPitch - 1;

ind = find(layer==1);
x(ind) = x(ind)+ crystalDepth;

phi = phi + double(rSector)*sectorPitch*pi/180;

% ---- Rotation by angle phi

xn = x;

xn = x.*cos(phi)-yo.*sin(phi);
yn = yo.*cos(phi)+x.*sin(phi);

xun = x.*cos(phi)-yu.*sin(phi);
yun = yu.*cos(phi)+x.*sin(phi);

xln = x.*cos(phi)-yl.*sin(phi);
yln = yl.*cos(phi)+x.*sin(phi);

end
