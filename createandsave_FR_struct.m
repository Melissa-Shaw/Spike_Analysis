function [Anaes_expFR] = createandsave_FR_struct(SA)

% set parameters
base_cond = 1;
sal_cond = 2;
tcb_cond = 3;

% create FR struct for anaesthetised recordings
Anaes_expFR.post = [];
Anaes_expFR.pre = [];
Anaes_expFR.TCB = [];
Anaes_expFR.expnum = [];
for e = 1:numel(SA)
  num_units = size(SA(e).M_cond_neuronFR,1);
  cond_idx = [zeros(num_units,1); ones(num_units,1)];
  exp_idx = repelem(SA(e).exp,num_units*2); % *2 to account for repeat pre values twice for compare with control and tcb
  
  Anaes_expFR.TCB = [Anaes_expFR.TCB; cond_idx];
  Anaes_expFR.post = [Anaes_expFR.post; SA(e).M_cond_neuronFR(:,sal_cond); SA(e).M_cond_neuronFR(:,tcb_cond)];
  Anaes_expFR.pre = [Anaes_expFR.pre; SA(e).M_cond_neuronFR(:,base_cond); SA(e).M_cond_neuronFR(:,base_cond)]; % repeat pre values twice for compare with control and tcb
  Anaes_expFR.expnum = [Anaes_expFR.expnum; exp_idx'];
end

Anaes_expFR.TCB = logical(Anaes_expFR.TCB);

% save firing rate matfile
save('E:\ms1121\Analysis Testing\Anaes_expFR.mat','Anaes_expFR');

end