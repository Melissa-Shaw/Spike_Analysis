
function [neuronFRcond,PopFRchangeP,FRchangeP,saline_spikerate,tcb_spikerate,minspikerate] = single_exp_spike_analysis(db, exp, spikestruct, saline_cond, tcb_cond, num_cond, sal_sect, tcb_sect, make_plots)

% extract relevant data
baseline_spikerate = spikestruct.spikerate(:,baseline_cond);
saline_spikerate = spikestruct.spikerate(:,saline_cond);
tcb_spikerate = spikestruct.spikerate(:,tcb_cond);
minspikerate = baseline_spikerate >= 0.05 & saline_spikerate >= 0.05 & tcb_spikerate >= 0.05; % cells with minimum spk rate of 0.05 spk/s in pre baseline and 0.05 spk/s in post baseline
cond_pop_rate = spikestruct.condpopspikevector;
num_units = numel(spikestruct.clusteridx);
spike_raster = spikestruct.raster;
timepoints = db(exp).timepointsms;
timepoints(1) = 1; % change from 0 to 1 for indexing

% find spikes in bins of 30sec
spikes = cell(1,num_cond);
cond_spike_raster = cell(1,num_cond);
neuronFRcond = cell(num_units,num_cond);
for c = 1:num_cond
  spikes{c} = calc_running_sum(cond_pop_rate{c},30000); % 30000ms = 30s
  cond_spike_raster{c} = spike_raster(:,timepoints(c):timepoints(c+1));
  for n = 1:num_units
    neuronFRcond{n,c} = calc_running_sum(cond_spike_raster{c}(n,:),1000); % bins of 1 second (1000ms) for spikes/s
  end
end

% ranksum to compare between conditions
[PopFRchangeP,~,~] = signrank(spikes{saline_cond}(sal_sect(1)/30:sal_sect(2)/30), spikes{tcb_cond}(tcb_sect(1)/30:tcb_sect(2)/30));
FRchangeP = NaN(num_units,1);
for n = 1:num_units
    [FRchangeP(n),~,~] = signrank(neuronFRcond{n,saline_cond}(sal_sect(1):sal_sect(2)), neuronFRcond{n,tcb_cond}(tcb_sect(1):tcb_sect(2)));
end

if make_plots == true
  % add condition lines to drift plot
  uiopen([db(exp).dir '\Drift_plot_all_spikes.fig'], 1)
  xlabel('Time (s)')
  y1 = ylim;
  for c = 1:num_cond
    xline(db(exp).timepointsms(c)/1000,'r-','linewidth',2);
    text([db(exp).timepointsms(c)/1000], y1(2), db(exp).injection{c});
  end
  title([db(exp).animal ' ' db(exp).date], 'Interpreter', 'none')

  % Plot firing rate and fit regression lines to conditions
  start = 1; % first value
  f1 = figure;
  tiledlayout('flow')
  nexttile;
  hold on
  for c = 1:num_cond
    min_idx = [start:(numel(spikes{c})+start-1)]./2; % gives time idx values in minutes
    firing_rate = spikes{c}./30; % firing rate is spikes per second
    regressionplot(min_idx,firing_rate,f1);
    start = (max(min_idx)+0.5)*2;
  end
  hold off
  yl = ylim;
  for c = 1:num_cond
    xline(db(exp).timepointsms(c)/60000,'r-','linewidth',2);
    text([db(exp).timepointsms(c)/60000], y1(2), db(exp).injection{c});
  end
  title([db(exp).animal ' ' db(exp).date], 'Interpreter', 'none')
  title(['Firing rate bins exp ' num2str(exp) ' Dose = ' num2str(db(exp).dose) ' N = ' num2str(numel(spikestruct.clusteridx)) ' '  sprintf(' p = %.2f ', PopFRchangeP(:))])
  ylabel('Firing rate (spikes/s)')
  xlabel('Time (minutes)')
  set(gca, 'fontsize', 14)

  % plot absolute change in mean firing rate vs Pre (> 0.5 sp/s)
  nexttile
  change_spikerate = tcb_spikerate - saline_spikerate; % Change in mean firing rate
  semilogx(saline_spikerate(minspikerate), change_spikerate(minspikerate), 'ko')
  yline(0,'k--');
  xlabel('Spikes per second (Saline)')
  ylabel('Change in mean firing rate Post-Pre (sp/s)')

  % Plot percentage change FR vs pre (> 0.5 sp/s)
  nexttile
  change_spikerate = tcb_spikerate - saline_spikerate; % Change in mean firing rate
  percchange_spikerate = (change_spikerate./saline_spikerate)*100; % Percentage Change in FR
  %percchange_spikerate = ((tcb_spikerate./saline_spikerate)*100);%-100; % percentage change in FR
  semilogx(saline_spikerate(minspikerate), percchange_spikerate(minspikerate), 'ko', 'Linewidth', 2)
  yline(0,'k--');
  xlim([10e-1 10e1])
  ylim([-100 350])
  xlabel('Saline FR (sp/s)')
  ylabel('Change in mean firing rate Post/Pre (%)')

  % Plot colour raster for FR change
  nexttile
  rasterwindow = 10000;
  raster = compressMatrix(spike_raster, 1, rasterwindow)*rasterwindow;
  baselinemeanFR = NaN(num_units,1);
  baselinesubtractFR = NaN(num_units,size(raster,2));
  zFR = baselinesubtractFR;
  for n = 1:num_units
      baselinemeanFR(n) = sum(spike_raster(n,1:300000))/300; % 5 min baseline (spk/s)
      baselinesubtractFR(n, :) = raster(n,:)/baselinemeanFR(n);
      zFR(n,:) = zscore(baselinesubtractFR(n, :));
  end
  [~, idx] = sort(baselinemeanFR);

  imagesc(zFR(idx,:), [-1 5])
  colorbar
  y1 = ylim;
  for c = 1:num_cond
    xline(timepoints(c)/rasterwindow,'r-','linewidth',2);
    text([db(exp).timepointsms(c)/rasterwindow], 1, db(exp).injection{c}, 'Color', 'red', 'FontSize', 14)
  end
  title('FR (10 sec bins) normalised to mean spike rate of first 5 minutes')
  ylabel('units sorted top to bottom by FR low to high')

  % Plot depth data
  if isfield(spikestruct, 'depth')
      % Pre FR vs depth
      nexttile
      plot(saline_spikerate, spikestruct.depth, 'ko')
      xlabel('Baseline FR (sp/s)')
      ylabel('Channel number')
      set(gca, 'XScale', 'log')

      % Percentage change FR vs depth
      nexttile
      plot(percchange_spikerate, spikestruct.depth, 'ko')
      xlabel('% Change FR (sp/s)')
      ylabel('Channel number')

      % Display progress
      disp(['Exp: ' num2str(exp) ' spike analyses complete'])
  else
      disp(['Exp: ' num2str(exp) ' spike analyses complete - no depth info'])
  end
  
end

end