function [ ospa, ospa_t ] = compute_ospa(est, truth, model, disp_flag)     
 % calculate and display average OSPA error

 ospa_t = zeros(1, truth.K);
 for k = 1:truth.K
     ospa_t(k) = ospa_dist( truth.X{k}, est.X{k}, model.ospa_cutoff, model.ospa_order);
 end
 ospa = nanmean(ospa_t);
 
 if strcmp(disp_flag, 'disp')
     fprintf(2,'\n  Average OSPA error: %d \n\n', ospa);
 end
     
end