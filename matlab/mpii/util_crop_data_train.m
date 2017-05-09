function annolist2 = util_crop_data_train(p, annolist, img_id, rectidxs)

saveDir = p.saveDir;
refHeight = p.refHeight;
deltaCrop = p.deltaCrop;
bCropIsolated = p.bCropIsolated;

if (~exist(saveDir, 'dir'))
    mkdir(saveDir);
end

annolist2 = [];
ncrops = 1;

delta_x = deltaCrop;
delta_y = deltaCrop;

numImgs = 0;

for imgidx = 1:length(annolist)
    
    fprintf('.');
    
    rect = annolist(imgidx).annorect;
    rect2 = annolist(imgidx).annorect;
    img = imread([p.imageDir annolist(imgidx).image.name]);
        
    for ridx = 1:length(rect)
        
        if (~isempty(setdiff(ridx,rectidxs{imgidx})))
            continue;
        end
        
        numImgs = numImgs + 1;
        
        if ~isfield(rect(ridx), 'annopoints') || isempty(rect(ridx).annopoints)
            continue;
        end
        
        pointsAll = [];
        points = rect(ridx).annopoints.point;
        for pid = 1:length(points)
            pp = [points(pid).x points(pid).y];
            pointsAll = [pointsAll; pp];
        end
        
        minX = min(pointsAll(:,1));
        maxX = max(pointsAll(:,1));
        minY = min(pointsAll(:,2));
        maxY = max(pointsAll(:,2));
        
        if (refHeight > 0)
            sc = rect(ridx).scale*200/refHeight;
            %% rescale image
            img_sc = imresize(img,1/sc,'bicubic');
        else
            sc = 1.0;
            img_sc = img;
            delta_x = round(deltaCrop*rect(ridx).scale);
            delta_y = round(deltaCrop*rect(ridx).scale);
        end
        
        posX1 = round(minX/sc);
        posX2 = round(maxX/sc);
        posY1 = round(minY/sc);
        posY2 = round(maxY/sc);
            
        %% crop image
        x1_new = round(max(1, posX1 - delta_x));
        x2_new = round(min(size(img_sc, 2), posX2 + delta_x));
        
        y1_new = max(1, posY1 - delta_y);
        y2_new = min(size(img_sc, 1), posY2 + delta_y);
        
        if (bCropIsolated && length(rect)>1)
            %% compute the closest annotated joint positions of other people
            points2All = [];
            for ridx2 = [1:ridx-1 ridx+1:length(rect2)]
                if (isfield(rect2(ridx2),'annopoints') && ~isempty(rect2(ridx2).annopoints) && ...
                    isfield(rect2(ridx2).annopoints,'point') && ~isempty(rect2(ridx2).annopoints.point))
                    points2 = rect2(ridx2).annopoints.point;
                    for pid = 1:length(points2)
                        pp = [points2(pid).x points2(pid).y];
                        points2All = [points2All; pp];
                    end
                end
            end
            if (~isempty(points2All))
                points2All = points2All./sc;
                d = points2All(:,1) - posX1; idx = find(d<0);
                [~,id] = max(d(idx)); posX1other = points2All(idx(id),1);
                d = points2All(:,2) - posY1; idx = find(d<0);
                [~,id] = max(d(idx)); posY1other = points2All(idx(id),2);
                d = posX2 - points2All(:,1); idx = find(d<0);
                [~,id] = max(d(idx)); posX2other = points2All(idx(id),1);
                d = posY2 - points2All(:,2); idx = find(d<0);
                [~,id] = max(d(idx)); posY2other = points2All(idx(id),2);
                
                if (refHeight > 0)
                    delta2 = refHeight/200*10;
                else
                    delta2 = rect(ridx).scale*10;
                end
                if (~isempty(posX1other))
                    x1_new = round(max(x1_new, posX1other+delta2));
                end
                if (~isempty(posX2other))
                    x2_new = round(min(x2_new, posX2other-delta2));
                end
            end
        end
        imgCrop = img_sc(y1_new:y2_new, x1_new:x2_new,1);
        imgCrop(:,:,2) = img_sc(y1_new:y2_new, x1_new:x2_new,2);
        imgCrop(:,:,3) = img_sc(y1_new:y2_new, x1_new:x2_new,3);
        
        %% save image
        fname = [saveDir '/im' padZeros(num2str(img_id(imgidx)),5) '_' num2str(ridx) '.png'];
        imwrite(imgCrop, fname);
        image_size = [size(imgCrop, 1), size(imgCrop, 2)];
        
        fnameT = [saveDir '/T_' padZeros(num2str(img_id(imgidx)),5) '_' num2str(ridx)];
        T = sc*[1 0 x1_new; 0 1 y1_new; 0 0 1];
        save(fnameT,'T');
        
        %% transfer annotations
        for pid = 1:length(points)
            points(pid).x = (points(pid).x)/sc - x1_new + 1;
            points(pid).y = (points(pid).y)/sc - y1_new + 1;
        end
        
        rect(ridx).annopoints.point = points;
        rect(ridx).x1 = rect(ridx).x1/sc - x1_new + 1;
        rect(ridx).y1 = rect(ridx).y1/sc - y1_new + 1;
        rect(ridx).x2 = rect(ridx).x2/sc - x1_new + 1;
        rect(ridx).y2 = rect(ridx).y2/sc - y1_new + 1;
        
        rect(ridx).objpos.x = rect(ridx).objpos.x/sc - x1_new + 1;
        rect(ridx).objpos.y = rect(ridx).objpos.y/sc - y1_new + 1;
        
        if (isempty(annolist2))
            annolist2.image.name = fname;
            annolist2.imgnum = 1;
            annolist2.annorect = rect(ridx);
            annolist2.image_size = image_size;
        else
            ncrops = ncrops + 1;
            annolist2(ncrops).image.name = fname;
            annolist2(ncrops).imgnum = ncrops;
            annolist2(ncrops).annorect = rect(ridx);
            annolist2(ncrops).image_size = image_size;
        end
        
    end
    
    if (~mod(imgidx, 100))
        fprintf(' %d/%d\n',imgidx,length(annolist));
    end
    
end
fprintf('\ndone\n');

%assert(numImgs == length(annolist2));

end
