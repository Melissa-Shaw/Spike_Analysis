%% running sum
  % inputs --> data vector and size of window to be summed

function [run_sum] = calc_running_sum(datafile,window)   
    trim_data = datafile(1:floor(numel(datafile)/window)*window); % trim to a length thats a multiple of the window required
    reshape_data = reshape(trim_data(:),1e3,[]);
    run_sum = sum(reshape_data);
end


    %run_sum = NaN(1,numel(trim_data)/window); % create NaN array of size
    %a = 1; % starting index point
    %for i=1:numel(run_sum)
    %  b = a+window-1; % end index point
    %  M = nansum(trim_data(a:b)); % sum for window
    %  run_sum(i) = M;
    %  a = b+1;
    %end
    %if isnan(run_sum)
    %  disp('NaN sum values present') % marking in case any windows were not calculated
    %end