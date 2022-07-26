function spikeanalysisPFC(db, experiments)
% Relies on db from makedb_TCB2
% First analysis script for PFC that creates figures for population:
% firing rate in 30s bins and fits regression lines to each condition
% change in firing rate between conditions of interest and tests for significance
% drift plot with conditions labelled
% baseline FR vs depth, percentage change to post drinking and depth


for exp = experiments
    topDir = 'R:\Neuropix\bd126\Analysis';
    method  = 'KS25';  % KS25 or manual
    expdir = [topDir '\Exp ' num2str(exp) ' ' db(exp).animal ' ' db(exp).injection{1:end-1}];
    baselineconditions = db(exp).cond;
    allconditions = 1:numel(db(exp).injection)-1;
    if strcmp(method, 'KS25')
        load([expdir '\spikestruct'])
    else
        load([expdir '\spikestruct'])
    end
    if isempty(spikestruct.clusteridx)
        disp(['Exp: ' num2str(exp) ' has no good units or is unsorted'])
        continue
    end
    %% save drift plot with lines in experiment directory
    if ~isfile([expdir '\Drift plot']) & ~strcmp(method, 'KS25')
        uiopen([db(exp).dir '\Drift_plot_all_spikes.fig'], 1)
        xlabel('Time (s)')
        y1 = ylim;
        for cond = 1:numel(db(exp).timepointsms)-1
            hold on, plot([1,1]*(db(exp).timepointsms(cond))/1000, [0, y1(2)], 'r-', 'linewidth', 2)
            text([db(exp).timepointsms(cond)/1000], y1(2), db(exp).injection{cond})
        end
        title([db(exp).animal ' ' db(exp).date], 'Interpreter', 'none')
        hgsave(gcf, [expdir '\Drift plot'])
        close
    end
    %% find unit parameters
    avspikerate = mean(spikestruct.spikerate, 2); % average spikerate across all conditions for each cell
    baselinespikerate = spikestruct.spikerate(:, baselineconditions(1));
    minspikerate = spikestruct.spikerate(:, baselineconditions(1)) >= 0.05 & spikestruct.spikerate(:, baselineconditions(2)) >= 0.05; % cells with minimum spk rate of 0.5 spk/s in pre baseline and 0.05 spk/s in post baseline
    %     isolationscore = spikestruct.unitquality(:,1)>= 20; % Isolation score >= 20
    %     refractoryperiod = spikestruct.unitquality(:,2)< 0.2; % contamination of 2ms refractory period relative to baseline < 0.20
    %     clusteridxgood = spikestruct.clusteridx(find(minspikerate & isolationscore & refractoryperiod)); % unit numbers of units meeting all criteria
    
    % find 30s windows of spikes per millisecond for each condition (population)
    windows = 3e4; % 30 second windows in ms samples
    for condition = allconditions % loop on conditions
        for nbins = 1:numel(db(exp).timepointsms(condition)+1:windows:db(exp).timepointsms(condition+1))-1 % loop on number of 30s bins in condition
            windowedges = db(exp).timepointsms(condition)+1:windows:db(exp).timepointsms(condition+1); % 30s windows
            % spikes1{condition} = sum(reshape(spikestruct.populationrate(floor(db(exp).timepointsms(condition)+1):floor(db(exp).timepointsms(condition+1)/3e4)*3e4), 3e4, [])',2);
            spikes{condition}(:,nbins) = sum(spikestruct.populationrate(floor(windowedges(nbins)):floor(windowedges(nbins+1)))); % number of spikes in each window
            for clu = 1:numel(spikestruct.clusteridx)
                neuronFRcond{clu, condition}(:,nbins) = sum(spikestruct.raster(clu, floor(windowedges(nbins)):floor(windowedges(nbins+1))),2);
            end
        end
    end
    
    % ranksum to compare between conditions
    [PopFRchangeP, h, stats] = ranksum(spikes{baselineconditions(1)}, spikes{baselineconditions(2)});
    for clu = 1:numel(spikestruct.clusteridx)
        [FRchangeP(clu,1),h,stats] = ranksum(neuronFRcond{clu,baselineconditions(1)}, neuronFRcond{clu,baselineconditions(2)});
    end
    save([expdir '\FRchange'], 'PopFRchangeP', 'spikes', 'FRchangeP', 'neuronFRcond', '-v7.3')
    
    % Plot firing rate and fit regression lines to conditions
    for condition = allconditions
        cumulativebinscond(condition) = sum(cellfun(@numel,spikes(1:condition))); % cumulative number of 30s bins for all conditions
    end
    
    f1 = figure;
    hold on
    for condition = allconditions
        for cumulativebins = sum(cellfun(@numel,spikes(1:condition)))
            if condition == 1
                regressionplot((1:cumulativebins)/2, spikes{condition}/30, f1)
            else
                regressionplot((cumulativebinscond(condition-1)+1:cumulativebinscond(condition))/2, spikes{condition}/30, f1) % spikes in 30s bins plotted in minutes(bins/2) and spikes per second(/30)
            end
        end
    end
    y1 = ylim;
    for condition = allconditions
        hold on,  plot([1,1]*numel(horzcat(spikes{1:condition-1}))/2, [0, y1(2)],'r-', 'linewidth', 2) % insert lines for conditions
        text(numel(horzcat(spikes{1:condition-1}))/2, y1(2), db(exp).injection{condition}, 'fontsize', 14) % insert condition labels
    end
    title(['Firing rate bins exp ' num2str(exp) ' Dose = ' num2str(db(exp).dose) ' N = ' num2str(numel(spikestruct.clusteridx)) ' '  sprintf(' p = %.2f ', PopFRchangeP(:))])
    ylabel('Firing rate (spikes/s)')
    xlabel('Time (minutes)')
    set(gca, 'fontsize', 14)
    hgsave(gcf, [expdir '\Firing rate bins goodunits - ' method])
    clear p
    close all
    
    % Plot spike rate between conditions (> 0.5 sp/s)
    
    %     [ratecorr, PvalRate] = corr(spikestruct.spikerate(minspikerate, baselineconditions(1)),...
    %         spikestruct.spikerate(minspikerate, baselineconditions(2)), 'rows', 'complete', 'type', 'Spearman');
    plot(spikestruct.spikerate(minspikerate, baselineconditions(1)), spikestruct.spikerate(minspikerate, baselineconditions(2)), 'o')
    %     title(sprintf('r = %.2f, p = %.3f', ratecorr, PvalRate))
    xlabel(['Spikes per second (' db(exp).injection{baselineconditions(1)} ')'])
    ylabel(['Spikes per second (' db(exp).injection{baselineconditions(2)} ')'])
    set(gca, 'XScale', 'log', 'YScale', 'log')
    hold on, plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], 'g--')
    hgsave(gcf, [expdir '\Spike rate ' db(exp).injection{baselineconditions(1)} ' vs ' db(exp).injection{baselineconditions(2)} ' exp ' num2str(exp) ' - ' method])
    close
    
    % plot absolute change in mean firing rate vs Pre (> 0.5 sp/s)
    SpikeRateChange = spikestruct.spikerate(:,baselineconditions(2)) - spikestruct.spikerate(:,baselineconditions(1)); % Change in mean firing rate
    %     [ratechangecorr, PvalRatechange] = corr(spikestruct.spikerate(minspikerate,baselineconditions(1)), SpikeRateChange(minspikerate),...
    %         'rows', 'complete', 'type', 'Spearman');
    semilogx(spikestruct.spikerate(minspikerate, baselineconditions(1)), SpikeRateChange(minspikerate), 'ko')
    %   text([spikestruct.spikerate(cond(1),spikestruct.clusteridx)],[SpikeRateChange], num2str(spikestruct.clusteridx1(1:end)))
    hold on, plot([min(xlim) max(xlim)], [0 0], 'k--')
    % title(['Exp ' num2str(exp) sprintf(' r = %.2f, p = %.3f', ratechangecorr(exp), PvalRatechange(exp))])
    xlabel(['Spikes per second (' db(exp).injection{baselineconditions(1)} ')'])
    ylabel('Change in mean firing rate Post-Pre (sp/s)')
    hgsave(gcf, [expdir '\Change in FR vs Pre Exp ' num2str(exp) ' - ' method])
    close
    
    % Plot percentage change FR vs pre (> 0.5 sp/s)
    SpikeRateChangepercent = ((spikestruct.spikerate(:, baselineconditions(2))./spikestruct.spikerate(:, baselineconditions(1)))*100)-100; % Percentage Change in FR
    %[ratechangepercentcorr, ~] = corr(spikestruct.spikerate(minspikerate, baselineconditions(1)), SpikeRateChangepercent(minspikerate),...
    % 'rows', 'complete', 'type', 'Spearman');
    semilogx(spikestruct.spikerate(minspikerate, baselineconditions(1)), SpikeRateChangepercent(minspikerate), 'ko', 'Linewidth', 2)
    xlim([10e-1 10e1])
    ylim([-100 350])
    hold on, plot([min(xlim) max(xlim)], [0 0], 'k--')
    % title(['Exp ' num2str(exp) sprintf(' r = %.2f, p = %.3f', ratechangepercentcorr, PvalRatechangepercent)])
    xlabel('Baseline FR (sp/s)')
    ylabel('Change in firing rate Post/Pre (%)')
    hgsave(gcf, [expdir '\Change in FR (%) vs Pre Exp ' num2str(exp) ' - ' method])
    close
    
    % plot colour raster for FR change
    rasterwindow = 10000;
    raster = compressMatrix(spikestruct.raster, 1, rasterwindow)*rasterwindow;
    for i = 1:size(spikestruct.raster, 1)
        baselinemeanFR(i) = sum(spikestruct.raster(i,1:300000))/300; % 5 min baseline (spk/s)
        baselinesubtractFR(i, :) = raster(i,:)/baselinemeanFR(i);
        zFR(i,:) = zscore(baselinesubtractFR(i, :));
    end
    [~, idx] = sort(baselinemeanFR);
    
    imagesc(zFR(idx,:), [-1 5])
    colorbar
    y1 = ylim;
    hold on
    for cond = 1:numel(db(exp).timepointsms)-1
        hold on, plot([1,1]*(db(exp).timepointsms(cond))/rasterwindow, [0, y1(2)], 'r-', 'linewidth', 2)
        text([db(exp).timepointsms(cond)/rasterwindow], 1, db(exp).injection{cond}, 'Color', 'red', 'FontSize', 14)
    end
    title('FR (10 sec bins) normalised to mean spike rate of first 5 minutes')
    ylabel('units sorted top to bottom by FR low to high')
    hgsave(gcf,[expdir '\FR change colourplot normalised to baseline - ' method])
    
    if isfield(spikestruct, 'depth')
        % Pre FR vs depth
        if strcmp(method, 'KS25')
            plot(spikestruct.spikerate(:,baselineconditions(1)), spikestruct.depth(spikestruct.goodunits), 'ko')
        else
            plot(spikestruct.spikerate(:,baselineconditions(1)), spikestruct.depth, 'ko')
        end
        xlabel('Baseline FR (sp/s)')
        ylabel('Channel number')
        set(gca, 'XScale', 'log')
        hgsave(gcf,[expdir '\FR vs Channel number - ' method])
        close
        
        % Percentage change FR vs depth
        if strcmp(method, 'KS25')
            plot(SpikeRateChangepercent, spikestruct.depth(spikestruct.goodunits), 'ko')
        else
            plot(SpikeRateChangepercent, spikestruct.depth, 'ko')
        end
        xlabel('% Change FR (sp/s)')
        ylabel('Channel number')
        hgsave(gcf,[expdir '\(%) Change FR vs Channel number - ' method])
        close
        disp(['Exp: ' num2str(exp) ' spike analyses complete'])
    else
        disp(['Exp: ' num2str(exp) ' spike analyses complete - no depth info'])
    end
    % clear FRchangeP neuronFRcond baselinesubtractFR
    clearvars -except topDir db
end % loop on db entries
end

