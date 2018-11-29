close all; clear all;
sortby="none";
i_common_code;
i_get_genes_residual_cv_value;
clearvars -except gl123 res_cv
load res_cv_neuron.mat T_res
[g,i,j]=intersect(gl123,T_res.Gene);
cv1=res_cv(i);
cv2=T_res.Residual_CV(j);
figure;
scatter(cv1,cv2)
corrmyown(cv1,cv2)

grid on
xlabel('Residual CV in LCL')
ylabel('Residual CV in neurons')
dt = datacursormode;
dt.UpdateFcn = {@i_myupdatefcn,g};


