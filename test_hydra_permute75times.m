function [] = test_hydra_permute75times
clear all;clc
num_perm = 75;

for ii = 1:num_perm,
%% NC (~50%) and realpatients (~50% of NC) randomly selected
[~, ARI] = run_hydra_experiment_csv_partialNCPatients('FILE.csv','.','C',0.25,'kmin',1,'kmax',8,'init',3,'cvfold',10);
ARI_val(:,ii) = ARI; clear ARI
end
save ARI_75permuteNCPatient ARI_val


for ii = 1:num_perm,
%% NC (~50%) and PseudoPatient (remaining 50%) randperm
[~, ARI] = run_hydra_experiment_csv_NC('FILE.csv','.','C',0.25,'kmin',1,'kmax',8,'init',3,'cvfold',10);
ARI_val(:,ii) = ARI; clear ARI
end
save ARI_75permuteNC_PseudoPatient ARI_val

end
  