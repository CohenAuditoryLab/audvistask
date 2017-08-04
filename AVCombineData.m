function AVCombineData
%% File selection and loading

% Select files
disp('Select both DDM files');
ddmfiles = uipickfiles('Type', {'*.mat', 'MAT-files'}, 'NumFiles', 2);
disp(ddmfiles(1));
disp(ddmfiles(2));
disp('Select all data files');
allfiles = uipickfiles('Type', {'*.mat', 'MAT-files'}, 'NumFiles', 2);
disp(allfiles(1));
disp(allfiles(2));

% Load files
ddm1 = load(ddmfiles{1});
all1 = load(allfiles{1});
ddm2 = load(ddmfiles{2});
all2 = load(allfiles{2});

%extract important variables
meta_data = all1.meta_data;
data1 = all2.meta_data;
ddm1 = ddm1.ddm_data;
ddm2 = ddm2.ddm_data;
all1 = all1.all_data; %all1.task_table;
all2 = all2.task_table;

%% Find like elements in both arrays, sorting by visual mode

%get unique visual values
v1 = table2cell(all1(:,2));
modes1 = unique(v1);
v2 = table2cell(all2(:,2));
modes2 = unique(v2);
allmodes = union(modes1, modes2);

%initialize data structures
ddm_data = [];
all_data = {};

%add sorted data to new structures
for i = 1:length(allmodes)
    inds1 = find(strcmp(v1, allmodes{i}));
    inds2 = find(strcmp(v2, allmodes{i}));
    ddm_data = [ddm_data; ddm1(inds1, :); ddm2(inds2, :)];
    all_data = [all_data; all1(inds1, :); all2(inds2, :)];
end

%% Save combined data structures

data_folder = '/Users/briannakarpowicz/Documents/Cohen Lab/Auditory-Visual Task/Data/';
save_filename = ['Combined_Data_' meta_data.subject '_' meta_data.date '_' data1.date];

meta_data = [meta_data, data1];

%save matlab data tables
save([data_folder save_filename '.mat'], 'all_data', 'meta_data');
save([data_folder 'DDM_' save_filename '.mat'], 'ddm_data');

end
    