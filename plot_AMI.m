

function [Awake_expFR,Anaes_expFR] = plot_AMI()

% read in firing rate mat file for awake recordings
load('E:\ms1121\Analysis Testing\AllexpFR.mat')
Awake_expFR = AllexpFR;
clear AllexpFR

% load in firing rate mat file for anaes recordings
load('E:\ms1121\Analysis Testing\Anaes_expFR.mat')

% find AMI
Awake_expFR.AMI = (Awake_expFR.pre)./(Awake_expFR.post);
Awake_expFR.AMI = abs(log(Awake_expFR.AMI));
Anaes_expFR.AMI = (Anaes_expFR.pre)./(Anaes_expFR.post);
Anaes_expFR.AMI = abs(log(Anaes_expFR.AMI));

% plot AMI for awake and anaes recordings
figure
tiledlayout('flow')
nexttile
boxplot([Anaes_expFR.AMI(~Anaes_expFR.TCB) Anaes_expFR.AMI(Anaes_expFR.TCB)],...
  'Color','k','Symbol',''); % plots boxplot of perc change
h = findobj(gca,'Tag','Box');
fill_color = [1 0 0; 0 0 0];
for j=1:length(h)
  patch(get(h(j),'XData'),get(h(j),'YData'),fill_color(j,:),'FaceAlpha',.25,'EdgeColor',fill_color(j,:),'LineWidth',2);
end
ylim([0 10]);
xticklabels({'Anaes con' 'Anaes tcb'});
awake_tcb = Awake_expFR.AMI(Awake_expFR.TCB);
awake_con = Awake_expFR.AMI(~Awake_expFR.TCB);
axis('square');
box off

nexttile
boxplot([awake_con(1:numel(awake_tcb)) awake_tcb],...
  'Color','k','Symbol',''); % plots boxplot of perc change
h = findobj(gca,'Tag','Box');
fill_color = [1 0 0; 0 0 0];
for j=1:length(h)
  patch(get(h(j),'XData'),get(h(j),'YData'),fill_color(j,:),'FaceAlpha',.25,'EdgeColor',fill_color(j,:),'LineWidth',2);
end
ylim([0 10]);
xticklabels({'Awake con' 'Awake tcb'})
axis('square');
box off


% probability density comparison
f1 = figure;
h = histogram(Awake_expFR.AMI(Awake_expFR.TCB),'normalization','pdf');
awake_tcb = [h.BinEdges(1:end-1)' h.Values'];
h = histogram(Awake_expFR.AMI(~Awake_expFR.TCB),'normalization','pdf');
awake_con = [h.BinEdges(1:end-1)' h.Values'];
h = histogram(Anaes_expFR.AMI(Anaes_expFR.TCB),'normalization','pdf');
anaes_tcb = [h.BinEdges(1:end-1)' h.Values'];
h = histogram(Anaes_expFR.AMI(~Anaes_expFR.TCB),'normalization','pdf');
anaes_con = [h.BinEdges(1:end-1)' h.Values'];

figure
tiledlayout('flow')
ax1 = nexttile;
plot(awake_con(:,1),awake_con(:,2),'k');
hold on
plot(awake_tcb(:,1),awake_tcb(:,2),'r');
hold off
xlabel('Absolute rate modulation');
ylabel('Probability Density');
legend({'Control','TCB-2'});
title('Awake');
box off

ax2 = nexttile;
plot(anaes_con(:,1),anaes_con(:,2),'k');
hold on
plot(anaes_tcb(:,1),anaes_tcb(:,2),'r');
hold off
xlabel('Absolute rate modulation');
ylabel('Probability Density');
legend({'Control','TCB-2'});
title('Anaesthetised');
box off

linkaxes([ax1 ax2],'xy');
xlim([0 3]);
ylim([0 2]);

end



