
[Tdata, TSdata, info] = load_open_ephys_data('G:\pavlovian_cholinergic_cellbase\HDB36\190513a\100_CH9.continuous');

load('G:\pavlovian_cholinergic_cellbase\HDB36\190513a\TT3_3.mat')

figure; 
hold on; 
plot(TSdata, Tdata); 
plot(TS, ones(1,numel(TS)));




cellid = 'HDB36_190513a_3.3';


dbstop if error;
global RESDIR;
cellbase = whichcb;
cbSwitcher(cellid);
[r, s, tetrodename] = cellid2tags(cellid);   %#ok<*ASGLU> % get tetrode number
matname = cellid2fnames(cellid,'cont',tetrodename);   % .mat filename for LFP
[pathname, filename, extension] = fileparts(matname);   % parse filename

filename = ['CSC' num2str(tetrodename)];
cscname = fullfile(pathname,[filename '.ncs']);   % filename for the Neuralynx CSC file


read_openephys('datadir', pathname, 'resdir', 'D:\_MATLAB_DATA\POOLED\rawTracePOOLED',...
    'TTspec', 3);

