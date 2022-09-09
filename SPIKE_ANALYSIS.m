%% SETUP
addpath('X:\cortical_dynamics\Shared\Code\matlib\stats');
addpath('X:\cortical_dynamics\User\ms1121\Code\General\');
run('makedb_TCB2_MS');
clear AwakeV1 Batch1PFC Batch2PFC Batch2V1 Batch3PFC

% set empty parameters for groupl plots
BASE.allFR = []; BASE.M_FR = []; BASE.SD_FR = [];
SAL.allFR = []; SAL.M_FR = []; SAL.SD_FR = [];
TCB.allFR = []; TCB.M_FR = []; TCB.SD_FR = [];
LOW_TCB.allFR = []; LOW_TCB.M_FR = []; LOW_TCB.SD_FR = [];

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
  BASE.allFR = [BASE.allFR; baseFR]; BASE.M_FR = [BASE.M_FR; nanmean(log10(baseFR(baseFR>=0.01)))]; BASE.SD_FR = [BASE.SD_FR nanstd(log10(baseFR(baseFR>=0.01)))];
  SAL.allFR = [SAL.allFR; salFR]; SAL.M_FR = [SAL.M_FR; nanmean(log10(salFR(salFR>=0.01)))]; SAL.SD_FR = [SAL.SD_FR nanstd(log10(salFR(salFR>=0.01)))];
  TCB.allFR = [TCB.allFR; tcbFR]; TCB.M_FR = [TCB.M_FR; nanmean(log10(tcbFR(tcbFR>=0.01)))]; TCB.SD_FR = [TCB.SD_FR nanstd(log10(tcbFR(tcbFR>=0.01)))];
  LOW_TCB.allFR = [LOW_TCB.allFR; low_tcbFR]; LOW_TCB.M_FR = [LOW_TCB.M_FR; nanmean(log10(low_tcbFR(low_tcbFR>=0.01)))]; LOW_TCB.SD_FR = [LOW_TCB.SD_FR nanstd(log10(low_tcbFR(low_tcbFR>=0.01)))];
  
  
  % save summary figure
  %savefig(['X:\cortical_dynamics\User\ms1121\Analysis Testing\Anaes_Spiking_Figures\Exp_Spiking_Summaries\Exp_' num2str(exp) '_Summary.fig']);
  %disp(['Exp:' num2str(exp) ' summary figure saved.']);
  
  e = e+1;
  
  end
end
clear e COI base_cond baseFR low_tcbFR sal_cond salFR t tcb_cond tcb_low_cond tcbFR

%% GROUP PLOTS

% remove zeros for log firing rate analysis
BASE.allFR(BASE.allFR < 0.01) = NaN;
TCB.allFR(TCB.allFR < 0.01) = NaN;
SAL.allFR(SAL.allFR < 0.01) = NaN;
LOW_TCB.allFR(LOW_TCB.allFR < 0.01) = NaN;


% set figure
figure
tiledlayout(2,3)
set(gcf,'Color','w');
set(gcf,'Position',get(0,'Screensize'));

% plot change in FR over all neurons from all recordings
nexttile(1)
plot_scatter_changeFR_all(BASE.allFR,SAL.allFR,'k.');
[~,p] = ttest(log10(BASE.allFR),log10(SAL.allFR)); p = round(p,3);
xlabel('Baseline FR'); ylabel('Control FR'); title(['N: ' num2str(sum(~isnan(BASE.allFR) & ~isnan(SAL.allFR))) ' p = ' num2str(p)]);
nexttile(2)
plot_scatter_changeFR_all(BASE.allFR,TCB.allFR,'r.');
[~,p] = ttest(log10(BASE.allFR),log10(TCB.allFR)); p = round(p,3);
xlabel('Baseline FR'); ylabel('TCB-2 FR'); title(['N: ' num2str(sum(~isnan(BASE.allFR) & ~isnan(TCB.allFR))) ' p = ' num2str(p)]);
nexttile(3)
plot_scatter_changeFR_all(BASE.allFR,LOW_TCB.allFR,'b.');
[~,p] = ttest(log10(BASE.allFR),log10(LOW_TCB.allFR)); p = round(p,3);
xlabel('Baseline FR'); ylabel('Low TCB-2 FR'); title(['N: ' num2str(sum(~isnan(BASE.allFR) & ~isnan(LOW_TCB.allFR))) ' p = ' num2str(p)]);

% plot cumulative distributions
nexttile(4)
plot_cumulative_dist(BASE.allFR,'k','--');
hold on
plot_cumulative_dist(SAL.allFR,'k','-');
plot_cumulative_dist(LOW_TCB.allFR,'b','-');
plot_cumulative_dist(TCB.allFR,'r','-');
hold off
legend({'Baseline' 'Control' 'Low TCB-2' 'TCB-2'},'Location','northeast');

% plot log gamma distributions
nexttile(5)
plot_gamma_dist(BASE.allFR,'k--');
hold on
plot_gamma_dist(SAL.allFR,'k');
plot_gamma_dist(LOW_TCB.allFR,'b');
plot_gamma_dist(TCB.allFR,'r');
hold off
legend({'Baseline' 'Control' 'Low TCB-2' 'TCB-2'},'Location','northeast');

% plot mean FR preVpost
%nexttile(7)
%plot_scatter(BASE.M_FR,SAL.M_FR,'k.');
%xlim([-2 2]); ylim([-2 2]);
%hold on
%plot_scatter(BASE.M_FR,TCB.M_FR,'r.');
%plot_scatter(BASE.M_FR,LOW_TCB.M_FR,'b.');
%xlabel('Mean logFR (before)'); ylabel('Mean logFR (after)');
%[~,p1] = ttest(BASE.M_FR,SAL.M_FR); p1 = round(p1,3);
%[~,p2] = ttest(BASE.M_FR,TCB.M_FR); p2 = round(p2,3);
%[~,p3] = ttest(BASE.M_FR,LOW_TCB.M_FR); p3 = round(p3,3);
%title({['BaseVsal: p = ' num2str(p1) ' BaseVtcb: p = ' num2str(p2)], ['BaseVlowtcb: p = ' num2str(p3)]});

% plot std FR preVpost
nexttile(6)
plot_scatter(BASE.SD_FR,SAL.SD_FR,'k.');
xlim([0 1.5]); ylim([0 1.5]);
hold on
plot_scatter(BASE.SD_FR,TCB.SD_FR,'r.');
plot_scatter(BASE.SD_FR,LOW_TCB.SD_FR,'b.');
xlabel('Std logFR (before)'); ylabel('Std logFR (after)');
[~,p1] = ttest(BASE.SD_FR,SAL.SD_FR); p1 = round(p1,3);
[~,p2] = ttest(BASE.SD_FR,TCB.SD_FR); p2 = round(p2,3);
[~,p3] = ttest(BASE.SD_FR,LOW_TCB.SD_FR); p3 = round(p3,3);
%title({['BaseVsal: p = ' num2str(p1) ' BaseVtcb: p = ' num2str(p2)], ['BaseVlowtcb: p = ' num2str(p3)]});


% LOCAL FUNCTIONS

function plot_scatter_changeFR_all(neuronFR_cond1,neuronFR_cond2,marker_style)
loglog(neuronFR_cond1,neuronFR_cond2,marker_style);
hline = refline(1,0);
hline.Color = [0.3 0.3 0.3];
hline.LineStyle = '--';
xlim([10^-2 10^2]);
ylim([10^-2 10^2]);
box off
axis square
end

function plot_cumulative_dist(FR,marker_colour,marker_style)
f = cdfplot(FR);
f.Color = marker_colour; f.LineStyle = marker_style;
xlabel('Firing Rate'); ylabel('Cumulative Probability'); title('');
set(gca,'Xscale','log');
xlim([10^-2 10^2]);
grid off
box off
axis square
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 1.0;
end
end

function plot_gamma_dist(FR,marker_style)
e = [-3:0.02:2];
FR_params = fitdist(FR, 'Gamma');
plot(e, log(10)*logammaPDF(FR_params.b, FR_params.a,e*log(10)),marker_style,  'LineWidth', 1); % logammaPDF gives natural log, need log(10)*logammaPDF) to convert to log10
xlabel('Log Firing Rate')
ylabel('Probability')
xlim([-2 2])
box off
axis square
end

function plot_scatter(xdata,ydata,marker_style)
scatter(xdata,ydata,200,marker_style);
hline = refline(1,0);
hline.Color = [0.3 0.3 0.3];
hline.LineStyle = '--';
box off
axis square
end
