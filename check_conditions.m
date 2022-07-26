function [sal_cond] = check_conditions(animal,date)

% check for saline condition
saline_recordings = {'M210316_MS','160321';...
                     'M210319_MS','190321';...
                     'M210521_MS','210521';...
                     'M210610_MS','100621';...
                     'M210611_MS','110621'};         
sal_cond = 0;
for R = 1:size(saline_recordings,1)
  if strcmp(animal,saline_recordings(R,1)) && strcmp(date,saline_recordings(R,2))
    sal_cond = 2;
  end
end

end