function [SA,align_point] = align_and_trim_FR_EDIT(SA)

  start_drinking = NaN(1,numel(SA));
  length = NaN(1,numel(SA));

  % Trim start to align according to first tcb2 injection (low or high)
  for i=1:numel(SA)
    [~,tcb_low_cond] = check_conditions(SA.animal,SA.date,SA.num_cond);
    if tcb_low_cond > 0
      start_drinking(i) = SA(i).timepoints(tcb_low_cond);
    else
      start_drinking(i) = SA(i).timepoints(SA.cond(2));
    end
  end

  align_point = min(start_drinking);
  for i=1:numel(SA)
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