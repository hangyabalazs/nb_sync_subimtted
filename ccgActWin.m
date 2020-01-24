function ccgActWin
% CCG activation window calculator

load('D:\_MATLAB_DATA\POOLED\ccgPOOLED\CCG_matrices_POOLED.mat');
lagZ = round(size(CCR,2)/2);
time = -lagZ:1:lagZ;
sigIntervals = CCR>UCCR;
allInt = zeros(size(sigIntervals));
for iS=1:size(sigIntervals,1)
    currVal = [0 sigIntervals(iS,:) 0];
    [peakInt, s1, s2] = unifying_and_short_killer(currVal, 10, 10);
    if ~isempty(s1)
        if any(diff(sort([s1 s2]))>1)
            allInt(iS,:) = peakInt;
            
            [~,sInx] = min(abs(s1-lagZ));
            [~,eInx] = min(abs(s2-lagZ));
            
            start1 = s1(sInx)-lagZ;
            end1 = s2(eInx)-lagZ;
            if (end1-start1)>1
                startPos(iS) = start1;
                endPos(iS) = end1;
                actLength(iS) = endPos(iS)-startPos(iS);
                absLength1(iS) = abs(start1);
                absLength2(iS) = abs(end1);
                [groupID, Label] = groupCellPair(group1{iS}, group2{iS});
                groups(iS) = groupID;
                
                cc1 = regexprep(PairOfCells{iS, 1},'_',' ');
                cc2 = regexprep(PairOfCells{iS, 2},'_',' ');
                H1 = figure;
                hold on;
                plot(time(2:end-1), CCR(iS,:), 'Color', [0 0 0.8])
                plot(time, currVal, 'Color', [0.8 0 0]);
                plot(time(2:end-1), UCCR(iS,:), 'Color', [0.2 0.2 0.2])
                plot(time(2:end-1), peakInt, 'Color', [0 0.8 0]);
                y_lim = ylim;
                line([start1 start1], y_lim, 'Color', [0 0 1 0.5], 'LineWidth', 2);
                line([end1 end1], y_lim, 'Color', [0 0 1 0.5], 'LineWidth', 2);
                title([cc1 ' ' cc2 ' ' Label...
                    ' range= ' num2str(start1) ' - ' num2str(end1) ' ms']);
                xlim([-250 250])
                maximize_figure;
                saveas(H1, ['D:\_MATLAB_DATA\POOLED\ccgActPOOLED\ccgAct' num2str(iS) '.jpeg']);
                close(H1);
            end
        end
    end
end

actStats.cells = PairOfCells;
actStats.groupID = groups;
actStats.start = startPos;
actStats.end = endPos;
actStats.actLength = actLength;
actStats.burstLength = mean(actLength(groups==1 | groups==4));
actStats.burstSE = nanse(actLength(groups==1 | groups==4));

allAbs = [abs(startPos) abs(endPos)];
bAVG = mean(allAbs(groups==1 | groups==4));
seAVG = nanse(allAbs(groups==1 | groups==4));

keyboard;