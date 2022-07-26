function [LFP] = get_LFP(db,exp,spikestruct,iCh)

%% Extract information of interest

LFPdata = spikestruct.LFP;

% extract key info
LFP.exp = exp;
if isfield(spikestruct,'animal')
  LFP.animal = spikestruct.animal;
  LFP.date = spikestruct.date;
end
if isfield(spikestruct,'meta')
  LFP.animal = spikestruct.meta.animal;
  LFP.date = spikestruct.meta.date;
end
%LFP.chan = db(exp).lfp;
LFP.dose = db(exp).dose;
LFP.gain = db(exp).LFPgain;

LFP.cond_timepoints = round(spikestruct.timepoints/1000);
if LFP.cond_timepoints(1) < 1 % correct for rounding error
  LFP.cond_timepoints(1) = 1;
end

% Set timepoints
cond_points = [1 numel(spikestruct.timepoints)];
start_time = round(spikestruct.timepoints(cond_points(1)));
end_time = round(spikestruct.timepoints(cond_points(2)));
if start_time < 1
  start_time = 1;
end
if end_time > numel(LFPdata{1})
  end_time = end_time-1;
end

if numel(LFP.cond_timepoints)==4
  LFP.start_drinking = LFP.cond_timepoints(2);
elseif numel(LFP.cond_timepoints)==6
  LFP.start_drinking = LFP.cond_timepoints(3);
else
  disp('Different number of conditions in recording')
end


% Set parameters
beta_freq = [12 30]; % beta frequency bounds
gamma_freq = [55 90]; % gamma frequency bounds
gaus = gausswin(60)/sum(gausswin(60)); % guassian filter
    
LFP.chan = db(exp).lfp(iCh);
    
% Extract raw LFP
LFP.raw = LFPdata{iCh}(start_time:end_time);

% Find saturations
[low_cut_off,up_cut_off] = manual_saturation(db);

if ~isnan(low_cut_off(exp)) | ~isnan(up_cut_off) 
  disp('Using manual saturation marking');
  up_sat_idx = find(LFP.raw > up_cut_off(exp));
  low_sat_idx = find(LFP.raw < low_cut_off(exp));
  sat_idx = [low_sat_idx up_sat_idx];
  LFP.sat_idx = unique(sat_idx);
else
  [LFP.sat_idx] = find_saturations(LFP.raw,15);
end

% Filter signal to extract beta and gamma
beta_LFP = leaveFrequencies(LFP.raw, 1e3, beta_freq(1), beta_freq(2)); % bandpass-filter frequency of interest, sampling rate at 1e3 for 1kHz
gamma_LFP = leaveFrequencies(LFP.raw, 1e3, gamma_freq(1), gamma_freq(2)); % bandpass-filter frequency of interest, sampling rate at 1e3 for 1kHz

 % Remove saturations
[LFP.raw] = remove_saturations(LFP.raw,LFP.sat_idx);
[beta_LFP] = remove_saturations(beta_LFP,LFP.sat_idx);
[gamma_LFP] = remove_saturations(gamma_LFP,LFP.sat_idx);

% Normalise the LFP (z-score)
[LFP.raw] = standardise_data(LFP.raw);
[beta_LFP] = standardise_data(beta_LFP);
[gamma_LFP] = standardise_data(gamma_LFP);

LFP.beta_st = beta_LFP;
LFP.gamma_st = gamma_LFP;

% Calculate the power
beta_LFP = beta_LFP(1:floor(numel(beta_LFP)/1e3)*1e3); % trim to a length thats a multiple of 1000
beta_LFP = reshape(beta_LFP(:),1e3,[]); 
LFP.beta = sum(beta_LFP.^2); % sum the squared values at 1Hz resolution
LFP.beta = LFP.beta/nanmean(LFP.beta(1:LFP.start_drinking));

gamma_LFP = gamma_LFP(1:floor(numel(gamma_LFP)/1e3)*1e3); % trim to a length thats a multiple of 1000
gamma_LFP = reshape(gamma_LFP(:),1e3,[]); 
LFP.gamma = sum(gamma_LFP.^2); % sum the squared values at 1Hz resolution

% Normalise the power by pre-baseline
LFP.beta = LFP.beta/nanmean(LFP.beta(1:LFP.start_drinking));
LFP.gamma = LFP.gamma/nanmean(LFP.gamma(1:LFP.start_drinking));

% Smooth the data across each minute
LFP.beta = nanconv(LFP.beta,gaus'); % convolution with gaussian filter
LFP.gamma = nanconv(LFP.gamma,gaus');

% Create placeholder variables
LFP.aligned_cond_timepoints = [];

disp(['LFP Exp: ' num2str(exp)  ' Cond_timepoints: ' num2str(cond_points) ' Complete']) % progress report


end