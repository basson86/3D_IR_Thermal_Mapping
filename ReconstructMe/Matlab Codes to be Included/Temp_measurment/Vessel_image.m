% 2012/11/25, testing the blood vessel segmention algorithm from :(2007)Physiology-Based Face Recognition in the thermal infrared spectrum
function [ImV]= Vessel_image(ImT)

% parameters for iteration 
niter= 5;
kappa= 20;
lambda= 0.1;
option= 1;

%         niter  - number of iterations.
%         kappa  - conduction coefficient 20-100 ?
%         lambda - max value of .25 for stability
%         option - 1 Perona Malik diffusion equation No 1
%                  2 Perona Malik diffusion equation No 2



% Applying anisotropic diffusion processing
diff = anisodiff(ImT, niter, kappa, lambda, option);



%%
se = strel('disk',8);
ImV = imtophat(diff,se);
