function [SA] = create_SA(db,exp,spikestruct,COI) % COI = conditions of interest (e.g. [1 2 3])

  % set parameters
  section = 5*60; % 5 minute section for mean FR
  gap = 5*60; % 5 minute gap into condition

  % extract relevant data
  SA.exp = exp;
  SA.animal = spikestruct.meta.animal;
  SA.date = spikestruct.meta.date;
  SA.region = db(exp).location;
  SA.num_units = numel(spikestruct.clusteridx);
  SA.clusteridx = spikestruct.clusteridx;
  SA.timepoints = round(db(exp).timepointsms/1000); 
  SA.timepoints(1) = 1; % change from 0 to 1 for indexing
  SA.num_cond = numel(db(exp).timepointsms)-1;
  SA.raster = spikestruct.raster;
  SA.cond = db(exp).cond;
  SA.design = spikestruct.design;
  
  % get neuron FR and population FR at 1Hz
  SA.neuronFR = NaN(SA.num_units,SA.timepoints(end));
  for n = 1:SA.num_units
    SA.neuronFR(n,:) = calc_running_sum(spikestruct.raster(n,:),1000); % bin of 1000ms for sp/s
  end
  SA.popFR = sum(SA.neuronFR);

  % get mean spikerate for each condition for each neuron
  SA.cond_neuronFR = cell(SA.num_units,SA.num_cond);
  SA.M_cond_neuronFR = NaN(SA.num_units,SA.num_cond);
  for c = 1:SA.num_cond
    for n = 1:SA.num_units
      SA.cond_neuronFR{n,c} = SA.neuronFR(n,SA.timepoints(c)+gap:SA.timepoints(c)+gap+section-1); % FR over section for indv neuron starting from gap after condition start
      SA.M_cond_neuronFR(n,c) = nanmean(SA.neuronFR(n,SA.timepoints(c)+gap:SA.timepoints(c)+gap+section-1)); % mean FR over section for indv neuron starting from gap after condition start
    end
  end
    
  % use step regression to find changepoint in FR
  SA.changepoint = NaN(SA.num_units,1);
  if isfield(spikestruct,'frameTimes')
      spont_neuronFR = [];
      SA.spont_timepoints = [1];
      for c = 1:SA.num_cond
          if isempty(spikestruct.frameTimes{c})
              spont_neuronFR = [spont_neuronFR SA.neuronFR(:,SA.timepoints(c):SA.timepoints(c+1)-1)]; % extract only spontaneous firing (exclude visual stimuli)
              SA.spont_timepoints = [SA.spont_timepoints size(spont_neuronFR,2)]; % extract new index timepoints for each extracted condition
          end
      end
  else
      spont_neuronFR = SA.neuronFR;
      SA.spont_timepoints = SA.timepoints;
  end
  
  for n = 1:SA.num_units
    [SA.changepoint(n)] = findchangepts(spont_neuronFR(n,:));
    if size(spont_neuronFR,2) < SA.changepoint(n)+gap+section-1 % if changepoint is too late for recording section to exist after it
      afterFR = nanmean(spont_neuronFR(n,SA.changepoint(n)+gap:end)); % just take changepoint to end of recording
    else
      afterFR = nanmean(spont_neuronFR(n,SA.changepoint(n)+gap:SA.changepoint(n)+gap+section-1)); % else take changepoint to end of recording section
    end
    if SA.changepoint(n)-gap-section < 1 % if changepoint is too early for recording section to exist before it
      beforeFR = nanmean(spont_neuronFR(n,1:SA.changepoint(n)-gap-1)); % just take from start to changepoint
    else
      beforeFR = nanmean(spont_neuronFR(n,SA.changepoint(n)-gap-section:SA.changepoint(n)-1)); % else take from recording section before changepoint
    end
    SA.changepoint_changeFR(n) = afterFR/beforeFR;
  end

end