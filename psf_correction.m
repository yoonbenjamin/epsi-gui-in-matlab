function [ img_abs_corr ] = psf_correction( img_abs )
% performs lucy-richardson correction
% % load('psf_map_011618.mat');
% % img_abs_corr=deconvlucy(img_abs,psf_map_011618);
   load('/Users/benjaminyoon/Desktop/PIGI folder/Projects/Project0 GUI/Code/Kidney Code/psf_map_022623.mat');
   img_abs_corr=deconvlucy(img_abs,urea_img_norm_xysym);
   psf_map = load('/Users/benjaminyoon/Desktop/PIGI folder/Projects/Project0 GUI/Code/Kidney Code/psf_map_022623.mat');
   psf_map_size = size(psf_map.urea_img_norm_xysym);
   disp(['PSF map size in MATLAB:', num2str(psf_map_size)]);
end