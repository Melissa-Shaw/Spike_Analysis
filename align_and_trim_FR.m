function [SA,align_point] = align_and_trim_FR(SA,tcb_cond)

  start_drinking = NaN(1,numel(SA));
  length = NaN(1,numel(SA));

  % Trim start to align according to start of drinking
  for i=1:numel(SA)
    start_drinking(i) = SA(i).timepoints(tcb_cond);
  end

  align_point = min(start_drinking);

  for i=1:numel(SA)
    %assert(numel(LFP(i).beta) == numel(LFP(i).gamma) && numel(LFP(i).beta) == size(LFP(i).specgram, 2))
    shift = start_drinking(i)-align_point;
    SA(i).timepoints = SA(i).timepoints-shift;
    SA(i).neuronFR = SA(i).neuronFR(:,shift+1:end);
    
    length(i) = size(SA(i).neuronFR,2);
  end
  
  % Trim end length to shortest recording
  trim_point = min(length);
  
  for i=1:numel(SA)
    SA(i).neuronFR = SA(i).neuronFR(:,1:trim_point);
  end
  
end