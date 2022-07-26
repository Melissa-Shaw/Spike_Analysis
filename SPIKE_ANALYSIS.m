%% SETUP
addpath('X:\cortical_dynamics\User\ms1121\Code');
run('makedb_TCB2_MS');
num_exp = numel(AnaesPFC);

% set empty parameters for groupl plots
all_baseFR = [];
all_salFR = [];
all_tcbFR = [];
all_low_tcbFR = [];

e = 1;
for exp = AnaesPFC
  
  % exclusions
  if exp ~= 142 && exp ~= 145 % exp 142 excluded for lack of units, exp 145 excluded for animal death
    
  % load spikestruct
  [spikestruct] = load_spikestruct('X:',db,exp);
  
  % extract conditions of interest
  base_cond = db(exp).cond(1);
  if numel(db(exp).cond)<3
    tcb_cond = db(exp).cond(2);
    tcb_low_cond = 0;
  else
    tcb_cond = db(exp).cond(3);
    tcb_low_cond = db(exp).cond(2);
  end
  [sal_cond] = check_conditions(db(exp).animal,db(exp).date);
  if sal_cond == 0
    COI = [base_cond tcb_low_cond tcb_cond];
  else
    COI = [base_cond sal_cond tcb_cond];
  end
  
  % create SA struct
  [SA(e)] = create_SA(db,exp,spikestruct,COI);
  
  % plot neuron FR data for each recording
  figure
  t = tiledlayout('flow');
  [baseFR,salFR,tcbFR,low_tcbFR] = plot_single_exp(SA(e),base_cond,sal_cond,tcb_cond,tcb_low_cond);
  set(gcf,'color','w');
  title(t,[SA(e).exp ' ' SA(e).animal ' ' SA(e).region ' Num_units: ' num2str(numel(SA(e).clusteridx))],Interpreter = 'none');
  
  % prepare neuron FR for each condition for group plots of all neurons
  all_baseFR = [all_baseFR; baseFR];
  all_salFR = [all_salFR; salFR];
  all_tcbFR = [all_tcbFR; tcbFR];
  all_low_tcbFR = [all_low_tcbFR; low_tcbFR];
  
  % save summary figure
  savefig(['X:\cortical_dynamics\User\ms1121\Analysis Testing\Anaes_Spiking_Figures\Exp_Spiking_Summaries\Exp_' num2str(exp) '_Summary.fig']);
  disp(['Exp:' num2str(exp) ' summary figure saved.']);
  
  e = e+1;
  
  end
end


%% GROUP PLOTS

% set figure
figure
tiledlayout('flow')

% plot change in FR over all neurons from all recordings
axs = [];
[ax1] = plot_scatter_changeFR_all(all_baseFR,all_salFR,'k','BaseVSal');
axs = [axs ax1];
[ax1] = plot_scatter_changeFR_all(all_baseFR,all_tcbFR,'r','BaseVTCB2');
axs = [axs ax1];
[ax1] = plot_scatter_changeFR_all(all_baseFR,all_low_tcbFR,'b','BaseVlow_TCB2');
axs = [axs ax1];
linkaxes(axs,'xy');
set(gcf,'Color','w');

%% COMPARE WITH AWAKE

% create and save FR struct for comparison with awake recordings
%[Anaes_expFR] = createandsave_FR_struct(SA);

% plot absolute modulation index for awake and anaes comparison
%plot_AMI();
%set(gcf,'color','w');


%% LOCAL FUNCTIONS

function [ax1] = plot_scatter_changeFR_all(neuronFR_cond1,neuronFR_cond2,marker_colour,manual_title)

ax1 = nexttile;
loglog(neuronFR_cond1,neuronFR_cond2,[marker_colour '.']);
%xlim([10^-2 10^2]);
%ylim([10^-2 10^2]);
hline = refline(1,0);
hline.Color = [0.3 0.3 0.3];
hline.LineStyle = '--';
xlabel('Firing rate (spikes/s) before')
ylabel('Firing rate (spikes/s) after')
%[R,corrP] = corrcoef(neuronFR_cond1,neuronFR_cond2);
title([manual_title ' N: ' num2str(sum(~isnan(neuronFR_cond1)&~isnan(neuronFR_cond2)))],'Interpreter','none');
box off
axis square

end
