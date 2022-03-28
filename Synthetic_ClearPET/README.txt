Scripts:

Paralel_MC_White_Image.m - MonteCarlo simulation of uniform distribution. As output you get the white image.
Reconstrucion_backprojection.m - (Filtered) Backprojection reconstrucion.
Reconstrucion_MLEM_like.m - Reconsturction by proposed method (MLEM).
Simulation_synthetic_measurement.m - Simulation of measurement. Crystals configuration cooresponds to ClearPET. 
White_image_expression.m - Generetes white image for one of two types of crystal configuration by closed-form expression.

-----------------------------------------------------------------------------
Functions:

ClearPET_gen_2d_single.m - ClearPET geometry.
CompToCartesian.m - Input: ClearPETs sectors, moduls, crystals and layers; Output: Cartesian coordinates of crystals.
nema_nu4_five.m - Model of NEMA NU 4-2008. Intersection with five cylinrars. Output: Cartesian coordinates of distribution of NEMA phantom.
nema_nu4_holes.m - Model of NEMA NU 4-2008. Intersection with two cylindric holes. Output: Cartesian coordinates of distribution of NEMA phantom.
rot_trian_1d_Rij.m - A close-form expression of rotation of triangular approximation of PDF.
rot_trian2.m - Summing over all pairs of crystals.


-----------------------------------------------------------------------------
Synthetic data:

comp_256_type1_075.m -The white image for type1 config. from the expression (image size:256×256)
comp_256_type1_075.m -The white image for type2 config. from the expression (image size:256×256)
comp_MC_256_type1_075.m -The white image for type1 config. by MC simulation (image size:256×256)
comp_MC_256_type1_075.m -The white image for type2 config. by MC simulation (image size:256×256)

Measurement_type1_fivespikes_075_30min - Data from synthetic measurement (config: type1; NEMA interserction: 5 ciricles(spikes); Numeber of coincidences : 30 minute equivalent)
Measurement_type2_fivespikes_075_30min - Data from synthetic measurement (config: type2; NEMA interserction: 5 ciricles(spikes); Numeber of coincidences : 30 minute equivalent)
Measurement_type1_twoholes_075_30min - Data from synthetic measurement (config: type1; NEMA interserction: 2 holes; Numeber of coincidences : 30 minute equivalent)
Measurement_type2_twoholes_075_30min - Data from synthetic measurement (config: type2; NEMA interserction: 2 holes; Numeber of coincidences : 30 minute equivalent)
