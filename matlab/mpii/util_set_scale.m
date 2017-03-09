function annolist = util_set_scale(annolist,ref_height)

if (nargin < 2)
    % reference height in px
    ref_height = 200;
end

HEAD_HEIGHT_RATIO = 1/8;

for imgidx = 1:length(annolist)
    if (isfield(annolist(imgidx), 'annorect') && ~isempty(annolist(imgidx).annorect))
        rect = annolist(imgidx).annorect;
        for ridx = 1:length(rect)
            if (isfield(rect(ridx), 'annopoints') && ~isempty(rect(ridx).annopoints))
                headSize = util_get_head_size(rect(ridx));
                sc = ref_height*HEAD_HEIGHT_RATIO/headSize;
                assert(sc<100 && sc>0.01);
                rect(ridx).scale = 1/double(sc);
            end
        end
        annolist(imgidx).annorect = rect;
    end
end