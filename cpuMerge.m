function cpuMerge
%CPUMERGE merges TAN psth figures for hit and fa

dbstop if error;
cb = 'HDB';
global RESDIR;
global PATH;
fs = filesep;

resdir = [RESDIR cb fs 'psth' cb fs 'figures' PATH fs];

% Figure1: -250 - 750
uiopen([resdir 'raw_psth_FB_allC_FA_original.fig'],1);
H_FA = gcf;
uiopen([resdir 'raw_psth_FB_allC_HIT_original.fig'],1);
H_HIT = gcf;
H_merge = figure;
A_merge = axes;
hold on
ln = findobj(H_FA,'type','line');
lnn = copyobj(ln,A_merge);
set(lnn,'Color',[145/255 56/255 0])
pt = findobj(H_FA,'type','patch');
ptn = copyobj(pt,A_merge);
set(ptn,'FaceColor',[145/255 56/255 0])
ln = findobj(H_HIT,'type','line');
lnn = copyobj(ln,A_merge);
set(lnn,'Color',[0 190/255 204/255])
pt = findobj(H_HIT,'type','patch');
ptn = copyobj(pt,A_merge);
set(ptn,'FaceColor',[0 190/255 204/255])
fnm = [resdir 'psth_merge_long.fig'];
fnmJ = [resdir 'psth_merge_long.jpeg'];
x_lim = [-250 750];
ylim([-1.5 10])
xticks([x_lim(1) 0 x_lim(1)*-1 x_lim(2)]);
yticks([0 5 10])
xlim(x_lim)
xlabel('Time from feedback (ms)')
ylabel('Firing rate (Hz)')
legend({'n=5'})
legend('boxoff')
setmyplot_tamas;
axis square;
saveas(H_merge, fnm);
saveas(H_merge, fnmJ);
close([H_HIT H_FA H_merge])

% Figure2: -20 - 80
uiopen([resdir 'raw_psth_FB_allC_FA_original.fig'],1);
H_FA = gcf;
uiopen([resdir 'raw_psth_FB_allC_HIT_original.fig'],1);
H_HIT = gcf;
H_merge = figure;
A_merge = axes;
hold on
ln = findobj(H_FA,'type','line');
lnn = copyobj(ln,A_merge);
set(lnn,'Color',[145/255 56/255 0])
pt = findobj(H_FA,'type','patch');
ptn = copyobj(pt,A_merge);
set(ptn,'FaceColor',[145/255 56/255 0])
ln = findobj(H_HIT,'type','line');
lnn = copyobj(ln,A_merge);
set(lnn,'Color',[0 190/255 204/255])
pt = findobj(H_HIT,'type','patch');
ptn = copyobj(pt,A_merge);
set(ptn,'FaceColor',[0 190/255 204/255])
fnm = [resdir 'psth_merge_short.fig'];
fnmJ = [resdir 'psth_merge_short.jpeg'];
x_lim = [-20 80];
ylim([0 50])
xticks([0 40 80]);
yticks([0 25 50])
xlim(x_lim)
xlabel('Time from feedback (ms)')
ylabel('Firing rate (Hz)')
legend({'n=5'})
legend('boxoff')
setmyplot_tamas;
axis square;
saveas(H_merge, fnm);
saveas(H_merge, fnmJ);
close([H_HIT H_FA H_merge])