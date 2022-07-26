function [p] = boxplot_conVtcb(con_perc_change,tcb_perc_change)
  p = round(signrank(con_perc_change,tcb_perc_change),3);
  boxplot([con_perc_change tcb_perc_change],'Color','k','Symbol',''); % plots boxplot of perc change
  %ylabel(['\Delta ' freq ' Power (%)']);
  h = findobj(gca,'Tag','Box');
  fill_color = [1 0 0; 0 0 0];
  for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),fill_color(j,:),'FaceAlpha',.25,'EdgeColor',fill_color(j,:),'LineWidth',2);
  end
  hold on
  %if p < 0.05
   %sigstar({[1 2]},p);
  %end
  x_val = ones(numel(con_perc_change),1);
  plot(x_val,con_perc_change,'o','MarkerFaceColor','k');
  x_val = x_val+1;
  plot(x_val,tcb_perc_change,'o','MarkerFaceColor','r');
  hold off
  %set(gca,'XTickLabel',{'Control' 'TCB-2'});
  if p < 0.001
   title(['p < 0.001'])
  else
    title(['p = ' num2str(p)]);
  end
  box off
end