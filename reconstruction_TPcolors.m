function reconstruction_TPcolors(cellids)
%RECONSTRUCTION_TPcolors   Anatomical reconstruction of electrode position.
%   RECONSTRUCTION_TPcolors(CELLIDS) plots the postition of CELLIDS on atlas images
%   base on the coordinates stored in CellBase (reference atlas: Franklin &
%   Paxinos, 3rd edition).

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   22-Sept-2013

%   Edit log: BH 9/22/13, TL 2/21/19

global RESDIR;
global cCode;
cellbase = whichcb;
datadir = [RESDIR cellbase '\reconstruction' cellbase '\'];

% Load atlas images (Franklin & Paxinos, 3rd edition)
I.p0p86 = imread([datadir 'Bp0p86.tif']);
corners.p0p86 = [46 962; 40 569];
corner_coordinates.p0p86 = [-4.75 4.75; 0.5 6];
dedicated_handle.p0p86 = 21;
isactive.p0p86 = false;

I.p0p74 = imread([datadir 'Bp0p74.tif']);
corners.p0p74 = [47 964; 47 575];
corner_coordinates.p0p74 = [-4.75 4.75; 0.5 6];
dedicated_handle.p0p74 = 22;
isactive.p0p74 = false;

I.p0p62 = imread([datadir 'Bp0p62.tif']);
corners.p0p62 = [51 967; 49 577];
corner_coordinates.p0p62 = [-4.75 4.75; 0.5 6];
dedicated_handle.p0p62 = 23;
isactive.p0p62 = false;

I.p0p50 = imread([datadir 'Bp0p50.tif']);
corners.p0p50 = [53 970; 53 584];
corner_coordinates.p0p50 = [-4.75 4.75; 0.5 6];
dedicated_handle.p0p50 = 24;
isactive.p0p50 = false;

I.m0p34 = imread([datadir 'Bm0p34.tif']);
corners.m0p34 = [89 1701; 75 1008];
corner_coordinates.m0p34 = [-4.75 4.75; 0.5 6];
dedicated_handle.m0p34 = 1;
isactive.m0p34 = false;

I.m0p46 = imread([datadir 'Bm0p46.tif']);
corners.m0p46 = [91 1703; 87 1020];
corner_coordinates.m0p46 = [-4.75 4.75; 0.5 6];
dedicated_handle.m0p46 = 2;
isactive.m0p46 = false;

I.m0p58 = imread([datadir 'Bm0p58.tif']);
corners.m0p58 = [89 1701; 78 1011];
corner_coordinates.m0p58 = [-4.75 4.75; 0.5 6];
dedicated_handle.m0p58 = 3;
isactive.m0p58 = false;

I.m0p70 = imread([datadir 'Bm0p70.tif']);
corners.m0p70 = [97 1709; 81 1014];
corner_coordinates.m0p70 = [-4.75 4.75; 0.5 6];
dedicated_handle.m0p70 = 4;
isactive.m0p70 = false;

I.m0p82 = imread([datadir 'Bm0p82.tif']);
corners.m0p82 = [115 1727; 100 1033];
corner_coordinates.m0p82 = [-4.75 4.75; 0.5 6];
dedicated_handle.m0p82 = 5;
isactive.m0p82 = false;

I.m0p94 = imread([datadir 'Bm0p94.tif']);
corners.m0p94 = [102 1714; 77 1010];
corner_coordinates.m0p94 = [-4.75 4.75; 0.5 6];
dedicated_handle.m0p94 = 6;
isactive.m0p94 = false;

I.m1p06 = imread([datadir 'Bm1p06.tif']);
corners.m1p06 = [99 1711; 100 1033];
corner_coordinates.m1p06 = [-4.75 4.75; 0.5 6];
dedicated_handle.m1p06 = 7;
isactive.m1p06 = false;

I.m1p22 = imread([datadir 'Bm1p22.tif']);
corners.m1p22 = [93 1705; 83 1016];
corner_coordinates.m1p22 = [-4.75 4.75; 0.5 6];
dedicated_handle.m1p22 = 8;
isactive.m1p22 = false;

I.m1p34 = imread([datadir 'Bm1p34.tif']);
corners.m1p34 = [80 1692; 69 1002];
corner_coordinates.m1p34 = [-4.75 4.75; 0.5 6];
dedicated_handle.m1p34 = 9;
isactive.m1p34 = false;

I.m1p46 = imread([datadir 'Bm1p46.tif']);
corners.m1p46 = [84 1696; 86 1019];
corner_coordinates.m1p46 = [-4.75 4.75; 0.5 6];
dedicated_handle.m1p46 = 10;
isactive.m1p46 = false;

I.m1p58 = imread([datadir 'Bm1p58.tif']);
corners.m1p58 = [89 1701; 83 1016];
corner_coordinates.m1p58 = [-4.75 4.75; 0.5 6];
dedicated_handle.m1p58 = 11;
isactive.m1p58 = false;

% TonicPhasic Groups
[phasicB, poissonL, tonic] = groupsTP(cellids);
% Order of plotting them on top of each other
cellsOrd = [cellids(poissonL) cellids(phasicB) cellids(tonic)];
cellids = cellsOrd;
groupID = [ones(sum(poissonL),1)*2; ones(sum(phasicB),1)*1; ones(sum(tonic),1)*3];

% Color code
mouse = getvalue('RatID',cellids);  % animal ID
nms = length(unique(mouse));   % number of mice
mouse2 = nan(size(mouse));
for k = 1:nms
    mouse2(mouse==min(mouse)) = k;
    mouse(mouse==min(mouse)) = Inf;
end
clr = colormap(jet(nms));   % unique color for each mouse
close;
mrk = 'Sp^hdo<>VX*+';   % unique marker for each mouse
colorcode = false;   % control color code
markercode = true;   % control marker code

% Cell coordinates
AP = getvalue('APpos',cellids);
DV = getvalue('DVpos',cellids);
L = getvalue('Lpos',cellids);

% ACG file for group identity
global PATH
cellbase = whichcb;
acgDir = [RESDIR cellbase '\acg' cellbase PATH '\'];
acgPath = [acgDir 'ACG_matrices_' cellbase '.mat'];
TPcolors = cCode;

% Put the cells on the map
for iC = 1:length(cellids)   % loop through the cells
    cellid = cellids{iC};   % current cell
    
    % Select the atlas image based on the AP position
    if AP(iC) < 0  % minus or plus: position rel. to Bregma
        morp = 'm';
    else
        morp = 'p';
    end
    i1 = floor(abs(AP(iC)/1000));  % integer part
    i2 = 100*abs(AP(iC)/1000-fix(AP(iC)/1000));  % fractional part (100th)
    if i2 < 10
        atlastag = [morp num2str(i1) 'p0' num2str(i2)];  % tag
    else
        atlastag = [morp num2str(i1) 'p' num2str(i2)];
    end
        
    % Open/activate figure
    figure(dedicated_handle.(atlastag))
    if ~isactive.(atlastag)
        imshow(I.(atlastag))
        isactive.(atlastag) = true;
        hold on
    end
    
    % Interpolate position
%     mlp = (L(iC) - corner_coordinates.(atlastag)(1,1)) / ...
%         (corner_coordinates.(atlastag)(1,2) - corner_coordinates.(atlastag)(1,1));
%     xpos = corners.(atlastag)(1,1) + mlp * ...
%         (corners.(atlastag)(1,2) - corners.(atlastag)(1,1));
    xpos = interp1(corner_coordinates.(atlastag)(1,:),corners.(atlastag)(1,:),-L(iC)/1000);
    ypos = interp1(corner_coordinates.(atlastag)(2,:),corners.(atlastag)(2,:),DV(iC)/1000);
    
    % Plot reconstructed position
    if colorcode
        cclr = clr(mouse2(iC),:);   % unique color for each mouse
    else
        lc = 22;   % number of identified neurons
        cclr = [0 0.8*(iC>lc) 0.8*(iC<=lc)];   % color according to identified/putative
    end
    
    cclr = TPcolors(groupID, :);
    
    if markercode
        cmrk = mrk(mouse2(iC));   % unique marker for each mouse
    else
        cmrk = 'o';
    end
    plot(xpos,ypos,cmrk,'MarkerSize',12,'MarkerFaceColor',cclr(iC,:),...
        'MarkerEdgeColor','k')
end

% Save
h =  findobj('type','figure');
for i = 1:length(h)
    saveas(gcf, [datadir 'Reconstructed' num2str(i) '.fig']);
    saveas(gcf, [datadir 'Reconstructed' num2str(i) '.jpeg']);
    close;
end