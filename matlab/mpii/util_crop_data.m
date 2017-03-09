function annolist2 = util_crop_data(annolist, saveDir, refHeight, img_id, rectidxs, groups, deltaCrop, deltaCrop2, padStr, bObjposOffset)

if (nargin < 3)
    refHeight = 200;
end

if (nargin < 7)
    deltaCrop = 150*refHeight/200;
end

if (nargin < 8)
    deltaCrop2 = 150*refHeight/200;
end

if (nargin < 10)
    bObjposOffset = true;
end

if (~exist(saveDir, 'dir'))
    mkdir(saveDir);
end

annolist2 = [];
ncrops = 1;

delta_x = deltaCrop;
delta_y = deltaCrop2;

numImgs = 0;

for imgidx = 1:length(annolist)
    
    fprintf('.');
    
    rect = annolist(imgidx).annorect;
    img = imread(annolist(imgidx).image.name);
    
    idxsGroup = zeros(length(rect),1); 
    scale = zeros(length(rect),1);
    pos = zeros(length(rect),2);
    
    for ridx = 1:length(rect)
        if (isfield(rect(ridx), 'annopoints') && ~isempty(rect(ridx).annopoints))
            gidx = 1;
            for idx = 1:length(groups{imgidx})
                if (ismember(ridx,groups{imgidx}{idx}))
                    gidx = idx;
                end
            end
            idxsGroup(ridx) = gidx;
            scale(ridx) = rect(ridx).scale*200/refHeight;
            pos(ridx,:) = [rect(ridx).objpos.x rect(ridx).objpos.y];
        end
    end
    
    idxsGroupUnique = unique(sort(idxsGroup));
        
    for i = 1:length(idxsGroupUnique)
        
        gidx = idxsGroupUnique(i);
        
        ridxs = find(idxsGroup == gidx);
        
        if (~isempty(setdiff(ridxs,rectidxs{imgidx})))
            continue;
        end
        
        numImgs = numImgs + 1;
        meanScale = mean(scale(ridxs));
        
        meanX = mean(pos(ridxs,1));
        meanY = mean(pos(ridxs,2));
        if (bObjposOffset)
            minX = min(pos(ridxs,1));
            minY = min(pos(ridxs,2));
            maxX = max(pos(ridxs,1));
            maxY = max(pos(ridxs,2));
        else
            pointsAll = [];
            rectGroup = rect(ridxs);
            for ridx = 1:length(rectGroup)
                points = rectGroup(ridx).annopoints.point;
                for pid = 1:length(points)
                    pp = round([points(pid).x points(pid).y]);
                    pointsAll = [pointsAll; pp];
                end
            end
            minX = min(pointsAll(:,1));
            maxX = max(pointsAll(:,1));
            minY = min(pointsAll(:,2));
            maxY = max(pointsAll(:,2));
        end
        
        posX1 = round(minX/meanScale);
        posX2 = round(maxX/meanScale);
        posY1 = round(minY/meanScale);
        posY2 = round(maxY/meanScale);
        
        pointsAll = [];
        rectGroup = rect(ridxs);
        
        for ridx = 1:length(rectGroup)
            points = rectGroup(ridx).annopoints.point;
            for pid = 1:length(points)
                pp = round([points(pid).x points(pid).y]/meanScale);
                pointsAll = [pointsAll; pp];
            end
        end
       
        %% rescale image
        if (refHeight > 0)
            img_sc = imresize(img,1/meanScale,'bicubic', 'antialiasing',false);
        else
            img_sc = img;
            delta_x = round(deltaCrop*mean(scale));
            delta_y = round(deltaCrop2*mean(scale));
        end
        
        delta_x1 = delta_x;
        delta_x2 = delta_x;
        delta_y1 = delta_y;
        delta_y2 = delta_y;
        
        %% crop image
        x1_new = round(max(1, posX1 - delta_x1));
        x2_new = round(min(size(img_sc, 2), posX2 + delta_x2));
        
        y1_new = max(1, posY1 - delta_y1);
        y2_new = min(size(img_sc, 1), posY2 + delta_y2);
        
        imgCrop = img_sc(y1_new:y2_new, x1_new:x2_new,1);
        imgCrop(:,:,2) = img_sc(y1_new:y2_new, x1_new:x2_new,2);
        imgCrop(:,:,3) = img_sc(y1_new:y2_new, x1_new:x2_new,3);
        
        %% save image
        fname = [saveDir '/im' padZeros(num2str(img_id(imgidx)),5) '_' num2str(gidx) '.png'];
        imwrite(imgCrop, fname);
        
        T = meanScale*[1 0 x1_new; 0 1 y1_new; 0 0 1];
        save([saveDir '/T_' padZeros(num2str(img_id(imgidx)),5) '_' num2str(gidx)],'T');
        
        for ridx = 1:length(rectGroup)
            
            %% transfer annotations
            points = rectGroup(ridx).annopoints.point;
            for pid = 1:length(points)
                points(pid).x = (points(pid).x)/meanScale - x1_new + 1;
                points(pid).y = (points(pid).y)/meanScale - y1_new + 1;
            end
            
            rectGroup(ridx).annopoints.point = points;
            rectGroup(ridx).x1 = rectGroup(ridx).x1/meanScale - x1_new + 1;
            rectGroup(ridx).y1 = rectGroup(ridx).y1/meanScale - y1_new + 1;
            rectGroup(ridx).x2 = rectGroup(ridx).x2/meanScale - x1_new + 1;
            rectGroup(ridx).y2 = rectGroup(ridx).y2/meanScale - y1_new + 1;
            
            rectGroup(ridx).objpos.x = rectGroup(ridx).objpos.x/meanScale - x1_new + 1;
            rectGroup(ridx).objpos.y = rectGroup(ridx).objpos.y/meanScale - y1_new + 1;
            
        end
        
        if (isempty(annolist2))
            annolist2.image.name = fname;
            annolist2.imgnum = 1;
            annolist2.annorect = rectGroup;
        else
            ncrops = ncrops + 1;
            annolist2(ncrops).image.name = fname;
            annolist2(ncrops).imgnum = ncrops;
            annolist2(ncrops).annorect = rectGroup;
        end
        
    end
    
    if (~mod(imgidx, 100))
        fprintf(' %d/%d\n',imgidx,length(annolist));
    end
    
end
fprintf('\ndone\n');

assert(numImgs == length(annolist2));

end