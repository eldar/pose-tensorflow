function [point, ind] = util_get_annopoint_by_id(points, id)
point = [];
for i=1:length(points)
    if (points(i).id == id)
        point = points(i);
        ind = i;
        return;
    end
end
end