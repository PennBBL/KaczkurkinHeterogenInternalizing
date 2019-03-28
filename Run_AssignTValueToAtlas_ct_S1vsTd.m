%Assigns the significant t-values to the JLF atlas.

%Specify the full path to the atlas. Should be .nii format.
Atlas_Path = '/data/jux/BBL/projects/pncHeterogeneity/images/mniJLF_LabelsWithWM.nii';

%Load the .csv file with the significant t-values and index numbers.
SigRegions = csvread('/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_ct_S1vsTd.csv');

%Define a vector of index numbers for the significant regions.
SigRegionIndex = SigRegions(:,end);

%Define a vector of t values for the significant regions.
TValue = SigRegions(:,1,:);

%Specify the path of the output file.
ResultantFile = '/data/jux/BBL/projects/pncHeterogeneity/images/ct_signifRegions_S1vsTd.nii';

%Run the AssignTValueToAtlas function.
AssignTValueToAtlas(Atlas_Path, SigRegionIndex, TValue, ResultantFile)
