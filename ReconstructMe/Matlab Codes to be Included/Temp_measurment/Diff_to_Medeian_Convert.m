function [ImT_Diff]= Diff_to_Medeian_Convert(ImT);
  
% Compute the mean value of the temperature over the domain of interest
ImT_median= median(median(ImT));
ImT_Diff=ImT - ImT_median;