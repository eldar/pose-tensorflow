function headSize = util_get_head_size(rect)

SC_BIAS = 0.6; % 0.8*0.75
headSize = SC_BIAS*norm([rect.x2 rect.y2] - [rect.x1 rect.y1]);

end