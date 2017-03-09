function preprocess_single(datasetDir, datasetName, saveDir, imageDir)

if (nargin < 2)
    datasetName = 'mpii_human_pose_v1_u12_1';
end

if (nargin < 3)
    saveDir = fullfile(datasetDir, 'cropped');
end

if (nargin < 4)
    imageDir = fullfile(datasetDir, 'images');
end

p = struct();

%crop_data(1,16000,1,400,130,1,0,1,130,0) % MPII single train
p.bTrain = 1;
p.refHeight = 400;
p.deltaCrop = 130;
p.bSingle = 1;
p.bCropIsolated = 1;
p.bMulti = 0;
p.bObjposOffset = 1;

p.datasetDir = datasetDir;
p.datasetName = fullfile(p.datasetDir, datasetName);
p.saveDir = saveDir;
p.imageDir = [imageDir '/'];

load(p.datasetName, 'RELEASE');
p.dataset = RELEASE;

annolist1 = crop_data(p);

%crop_data(1,16000,1,400,65,0,0,0,65,0) % multiperson as single people, train
p.deltaCrop = 65;
p.bSingle = 0;
p.bCropIsolated = 0;

annolist2 = crop_data(p);

annolist = horzcat(annolist1, annolist2);

annolistFullName = [p.saveDir '/annolist-full-h' num2str(p.refHeight)];
save(annolistFullName, 'annolist');

prepare_training_data(annolist, p.saveDir);

end
