%Cortical Thickness
%% Setup
cd /Users/efstathioskondylis/Desktop/CingulateProject/Imaging
Subjects = {'Sbj1','Sbj2','Sbj7'};
Map = 'DKTatlas40';
Gyri = {'caudalanteriorcingulate', 'isthmuscingulate',...
    'posteriorcingulate', 'rostralanteriorcingulate'};

%% Run
[v,l,ct] = read_annotation(fullfile(cd,'stats','lh.aparc.a2009s.annot'));