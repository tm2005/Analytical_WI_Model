Scripts:

const_list - ClearPET geometry (constants)
From_bin_to_cartisian - Same as CompToCartesian, but it also takes intersection with one possiable axial valus. Variable: "loc"
Reconstrucion_backprojection.m - (Filtered) Backprojection reconstrucion.
Reconstrucion_MLEM_like.m - Reconsturction by proposed method (MLEM).
Reconstruct_all_intersections - Reconsturction by proposed method (MLEM). Takes one measurement and take all 48 intersections and preforms proposed method on each of them.
    Afterwards, all reconstructed images are save at (by default) /set1/imgN.png. (N say which intersection)
White_image_expression.m - Generetes white image for one of two types of crystal configuration by closed-form expression.

-----------------------------------------------------------------------------
Functions:

ClearPET_gen_2d_single.m - ClearPET geometry.
CompToCartesian - Data (Events) for ClearPET are reprezented by Sector, Module, Layer, Crystal and Mechanical angle. For example
    Sector(55,1)=1, Sector(55,2)=12, Module(55,1) = 2, Module(55,2) = 3, Layer(55,1)=1, Layer(55,1)=2, Crystal(55,1)=38, Crystal(55,2) = 27, angle(55)=26.8 deg.
    That means that 55th captured events was between sectors 1 and 12. The one that hit sector 1 also hit module 2, crystal number 38 in layer 1. 
    The one that hit sector 2 also hit module 3, crystal number 27 in layer 2. This event occured when mechanical angle was at 26.8 deg.
    This function takes Sectors, Modules, Layers, Crystals and Mechanical angles on one measurement and turns them into Cartesian coordinates.
rot_trian_1d_Rij.m - A close-form expression of rotation of triangular approximation of PDF.
rot_trian2.m - Summing over all pairs of crystals.


-----------------------------------------------------------------------------
Synthetic data:

allzn.mat - List of all possible z-coordinates that coorespond to crystals center. They coorespond to different axial intersection. There are 48 such intersections.
          - ((-54.05) + 2.3*k) mm , k=0,1,2,3,...,47

comp_256_type1_075.mat -The white image for type1 config. from the expression (image size:256×256)
comp_256_type1_075.mat -The white image for type2 config. from the expression (image size:256×256)
comp_MC_256_type1_075.mat -The white image for type1 config. by MC simulation (image size:256×256)
comp_MC_256_type1_075.mat -The white image for type2 config. by MC simulation (image size:256×256)

Measurement_type1_fivespikes_075_30min.mat - Data from synthetic measurement (config: type1; NEMA interserction: 5 ciricles(spikes); Numeber of coincidences : 30 minute equivalent)
Measurement_type2_fivespikes_075_30min.mat - Data from synthetic measurement (config: type2; NEMA interserction: 5 ciricles(spikes); Numeber of coincidences : 30 minute equivalent)
Measurement_type1_twoholes_075_30min.mat - Data from synthetic measurement (config: type1; NEMA interserction: 2 holes; Numeber of coincidences : 30 minute equivalent)
Measurement_type2_twoholes_075_30min.mat - Data from synthetic measurement (config: type2; NEMA interserction: 2 holes; Numeber of coincidences : 30 minute equivalent)
