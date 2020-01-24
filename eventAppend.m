function eventAppend


for iC = 1:numel(allChAT)
cbSwitcher(allChAT{iC})
prealignSpikes(allChAT{iC},'FUNdefineEventsEpochs',@defineEventsEpochs_trialstart,...
'filetype','event','ifsave',1,'ifappend',1)
end