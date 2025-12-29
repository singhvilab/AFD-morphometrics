function allData=eteetsWrap()

%This is a wrapper function to analyze NREs shape related to the work of
%Teets et al from the Singhvi lab.
% uigetdir() should point to the directory containing the tif stacks to
% analyze
% output:
%   allData: numImages-by-1 cell containing computed
%   results for image in a table form

%% Select image directory and initialize data
selpath = uigetdir;
if selpath == 0
    allData = [];
    return
end
cd(selpath);
pd = dir('*.tif');
nIm = numel(pd);
allData = cell(nIm, 1);

%% Compute metrics for each image
for i = 1:nIm

    [Icleantrcrop, NWsk, DW, DistT, FinalT] = eteetsNeurons(pd(i).name);
    allData{i,1} = FinalT;
    writetable(DistT, 'distances.xlsx', 'sheet', i);
    fname = pd(i).name;
    sname = [fname(1:end-4) '.mat'];
    save(sname, 'Icleantrcrop', 'NWsk', 'DW');
end

%% Export data to Excel file to location of choice
AT = vertcat(allData{:});
file = uiputfile('*.*', 'Save Data...', 'myData.xlsx');
writetable(AT, file);