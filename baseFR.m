function baseFR

load('D:\_MATLAB_DATA\POOLED\acgPOOLED\ACG_matrices_POOLED.mat')


for iC=1:numel(cellids)
   currCell = cellids{iC};
   cbSwitcher(currCell);
   
   FR = firingrate_analysis(currCell);
   allFR(iC) = FR;
    
end



keyboard;