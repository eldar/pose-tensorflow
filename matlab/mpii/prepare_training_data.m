function prepare_training_data(annolist, saveDir)

zero_indexed_joint_ids = true;

pidxs = [0 2 4 5 7 9 12 14 16 17 19 21 22 23];
num_joints = length(pidxs);

load('parts.mat');

mkdir_if_missing(saveDir);

num_images = length(annolist);
channels = 3; % three channel images

dataset = struct('image',{}, 'size',{}, 'joints', {});

for i = 1:num_images
    if mod(i, 100) == 0
        fprintf('processing image %d/%d \n', i, num_images);
    end
    filename = annolist(i).image.name;
    
    joints = zeros(num_joints, 3);
    num_people = length(annolist(i).annorect);
    all_joints = cell(1,1);
    for k = 1:num_people
        joint_list = get_anno_joints(annolist(i).annorect(k), pidxs, parts);
        
        n = 0;
        for j = 1:num_joints
            jnt = joint_list(j, :);
            if ~isnan(jnt(1))
                n = n + 1;
                joints(n, :) = [j jnt];
            end
        end
        joints = int32(joints(1:n, :));
        if zero_indexed_joint_ids
            joints(:,1) = joints(:,1) - 1;
        end
        all_joints{k} = joints;
    end
    
    entry = struct;
    entry.image = filename;
    entry.size = [channels, annolist(i).image_size];
    entry.joints = all_joints;
    dataset(i) = entry;
end

out_filename = fullfile(saveDir, 'dataset.mat');
fprintf('Generated dataset definition file: %s\n', out_filename);
save(out_filename, 'dataset');

end

