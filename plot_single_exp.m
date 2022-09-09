function [Baseline_FR,Saline_FR,Tcb_FR,Tcb_low_FR] = plot_single_exp(SA,base_cond,sal_cond,tcb_cond,tcb_low_cond)

Baseline_FR = SA.M_cond_neuronFR(:,base_cond);
Tcb_FR = SA.M_cond_neuronFR(:,tcb_cond);

if sal_cond > 0
  Saline_FR = SA.M_cond_neuronFR(:,sal_cond);
else
  Saline_FR = NaN(size(Baseline_FR,1),1);
end
if tcb_low_cond > 0
  Tcb_low_FR = SA.M_cond_neuronFR(:,tcb_low_cond);
else
  Tcb_low_FR = NaN(size(Baseline_FR,1),1);
end

% plot raster of spikes at second frequency
ax1 = nexttile;
plotSpikeRaster(SA.neuronFR>0,'PlotType','vertline');
xlabel('Time (s)');
ylabel('Unit');
for c = 2:numel(SA.timepoints)-1
  xline(SA.timepoints(c),'r');
  txt = text(SA.timepoints(c),-0.3,SA.design(c),Interpreter = 'none');
  set(txt,'Rotation',90);
end

% plot log scatter of baslineVsaline neuron FR
if sal_cond > 0
  nexttile
  [r,~] = plot_neuronFR_scatter(Baseline_FR,Saline_FR,'k.');
  title(['Baseline vs Saline (r = ' num2str(r) ')']);
end

% plot log scatter of baselineVtcb neuron FR
nexttile
[r,~] = plot_neuronFR_scatter(Baseline_FR,Tcb_FR,'r.');
title(['Baseline vs TCB2 (r = ' num2str(r) ')']);

% plot log scatter of baslineVtcb_low neuron FR
if tcb_low_cond > 0  
  nexttile
  [r,~] = plot_neuronFR_scatter(Baseline_FR,Tcb_low_FR,'b.');
  title(['Baseline vs Low TCB2 (r = ' num2str(r) ')']);
end

% plot population rate in 1s bins
ax2 = nexttile;
plot(smooth(SA.popFR,100));
for c = 2:numel(SA.timepoints)
  xline(SA.timepoints(c),'r');
end
xlabel('Time (s)');
ylabel('Population FR');
linkaxes([ax1 ax2],'x');

% boxplot sum of squares for salineVtcb
if sal_cond > 0
  nexttile
  con = Saline_FR./Baseline_FR;
  tcb = Tcb_FR./Baseline_FR;
  boxplot_conVtcb(con,tcb);
  ylabel('Ratio of firing rate to baseline');
  set(gca,'XTickLabel',{'Saline' 'TCB-2'});
end

% boxplot sum of squares for tcb_lowVtcb
if tcb_low_cond > 0
  nexttile
  tcb_low = Tcb_low_FR./Baseline_FR;
  tcb = Tcb_FR./Baseline_FR;
  boxplot_conVtcb(tcb_low,tcb);
  ylabel('Ratio of firing rate to baseline');
  set(gca,'XTickLabel',{'TCB-2 Low' 'TCB-2'});
end

% plot changeFR around changepoint
nexttile
semilogy(SA.changepoint,SA.changepoint_changeFR,'k.');
xlim([0 SA.spont_timepoints(end)]);
xline(SA.spont_timepoints(3),'r');
txt = text(SA.spont_timepoints(3)-70,1,'TCB2');
set(txt,'Rotation',90);
xline(SA.spont_timepoints(2),'r');
if tcb_low_cond > 0
    %xline(SA.spont_timepoints(2),'r');
    txt = text(SA.spont_timepoints(2)-70,1,'Low TCB-2');
    %set(txt,'Rotation',90);
else
    txt = text(SA.spont_timepoints(2)-70,1,'Control');
end
set(txt,'Rotation',90);
%for c = 2:numel(SA.timepoints)
%  xline(SA.timepoints(c),'r');
%end
yline(1,'--');
ylabel('\Delta Firing Rate');
xlabel('Time (s)');
ylim([-1000 1000]);
%title('Step Changepoints');

%% Local Functions

function [r,corrP] = plot_neuronFR_scatter(x_Data,y_Data,marker_style)
  loglog(x_Data,y_Data,marker_style);
  xlim([10^-2 10^2]);
  ylim([10^-2 10^2]);
  hline = refline(1,0);
  hline.Color = [0.3 0.3 0.3];
  hline.LineStyle = '--';
  xlabel('Firing rate (spikes/s) before')
  ylabel('Firing rate (spikes/s) after')
  [R,corrP] = corr(x_Data,y_Data);
  r = R; %R(1,2); % R is a square matrix of x and y vs x and y, need only the r value for x vs y
  box off
  axis square
end


end