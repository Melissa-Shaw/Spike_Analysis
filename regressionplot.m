% Inputs: x, y - vectors
%figH is the handle of the figure of the regression line superimposed on
%the data
% plotSignificantOnly (optional) will plot only significant regression output if entered
% The function plots the regression of y on x.
function regressionplot(x, y, figH)
x = double(torow(x)');
y = double(torow(y)');
x(1:end, 2) = ones;
[b,bint,r,rint,stats] = regress(y,x);

figure(figH)
plot(x(1:end,1),y,'.')
xx = [min(x(1:end,1)), max(x(1:end,1))];
yy = [xx(1)*b(1) + b(2), xx(2)*b(1) + b(2)];

if stats(3)<.05
  plot(xx, yy, 'r')
else
  plot(xx, yy, 'b')
end
end
%% Test code:
%figH = figure; hold on
%x = randn(1, 100);
%y = 3.99*x + randn(1, 100);
%regressionplot(x, y, figH)
%regressionplot(x+7, y+73, figH)
