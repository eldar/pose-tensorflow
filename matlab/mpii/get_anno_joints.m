function [ joints ] = get_anno_joints( rect, pidxs, parts )

num_joints = length(pidxs);
joints = NaN(num_joints, 2);

points = rect.annopoints.point;
for j = 1:num_joints
    pidx = pidxs(j);
    annopoint_idxs = parts(pidx+1).pos;
    assert(annopoint_idxs(1) == annopoint_idxs(2));
    pt = util_get_annopoint_by_id(points, annopoint_idxs(1));
    if (~isempty(pt))
        joints(j, :) = [pt.x pt.y];
    end
end

end

