function annolist = crop_data(p)

fprintf('crop_data()\n');

bTrain = p.bTrain;
refHeight = p.refHeight;
deltaCrop = p.deltaCrop;
bSingle = p.bSingle;
bCropIsolated = p.bCropIsolated;
bMulti = p.bMulti;
bObjposOffset = p.bObjposOffset;
saveDir = p.saveDir;
DATASET = p.dataset;
annolist = DATASET.annolist;

func_crop_data_test = @util_crop_data;
func_crop_data_train = @util_crop_data_train;

fprintf('bTrain: %d\n',bTrain);
fprintf('refHeight: %d\n',refHeight);
fprintf('deltaCrop: %d\n',deltaCrop);
fprintf('bSingle: %d\n',bSingle);
fprintf('bCropIsolated: %d\n',bCropIsolated);
fprintf('bMulti: %d\n',bMulti);
fprintf('bObjposOffset: %d\n',bObjposOffset);

if (bTrain)
    split = 'train-v15';
else
    split = 'test-v15';
end

%{
if (bSingle)
    mode = 'singlePerson';
    rectidxs = DATASET.single_person;
else
    mode = 'multPerson';
    for imgidx = 1:length(DATASET.single_person)
        DATASET.mult_person{imgidx} = [DATASET.mult_person{imgidx} DATASET.borderline_person{imgidx}'];
    end
    rectidxs = DATASET.mult_person;
end
%}


if (bSingle)
    mode = 'singlePerson';
    rectidxs = DATASET.single_person;
else
    mode = 'multPerson';
    rectidxs = cell(size(DATASET.single_person));
    for imgidx = 1:length(DATASET.single_person)
        single_person = DATASET.single_person{imgidx};
        single_person = reshape(single_person, [1, length(single_person)]);
        rectidxs{imgidx} = setdiff(1:length(annolist(imgidx).annorect), single_person);
    end
end

annolistFullName = [saveDir '/annolist-' mode '-h' num2str(refHeight) '.mat'];
if exist(annolistFullName, 'file') == 2
    load(annolistFullName, 'annolist');
    return;
end

imgidxs1 = find(cellfun(@isempty,rectidxs) == 0);
imgidxs2 = find(DATASET.img_train == bTrain);
imgidxs = intersect(imgidxs1,imgidxs2);

imgidxs_sub = imgidxs;

% set scale
annolistSubset = util_set_scale(annolist(imgidxs_sub),200);

% crop images
if (~bTrain) % test
    assert(false);
    annolist = func_crop_data_test(annolistSubset, saveDir, refHeight, imgidxs_sub, rectidxs(imgidxs_sub), DATASET.groups(imgidxs_sub), deltaCrop, deltaCrop, 'symmetric', bObjposOffset);
else % train
    if (bMulti)
        assert(false);
        annolist = func_crop_data_test(annolistSubset, saveDir, refHeight, imgidxs_sub, rectidxs(imgidxs_sub), DATASET.groups(imgidxs_sub), deltaCrop, deltaCrop, 'symmetric', false);
    else
        annolist = func_crop_data_train(p, annolistSubset, imgidxs_sub, rectidxs(imgidxs_sub));
    end
end

save(annolistFullName, 'annolist');
