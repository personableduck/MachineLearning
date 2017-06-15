function [ centroidFrom, centroidTo ] = MatchPoints( centroidOrig,distanceThreshold,minIntensities )
%For each point, find closest point, check distance
        %Take lower pix val to be from traveled, take higher pix val to be to traveled

centroidFrom = [];
centroidTo = [];
while(size(centroidOrig,1)>1)
    currPointX = centroidOrig(1,1);
    currPointY = centroidOrig(1,2);
    otherX = centroidOrig(2:end,1);
    otherY = centroidOrig(2:end,2);
    [currIdx, currDist] = find_min_dist(currPointX,currPointY,otherX,otherY);
    currIdx = currIdx + 1; %The offset comparing otherX,otherY arrays to centroidOrig matrix
    
    %check distance
    if(currDist>distanceThreshold)
        centroidOrig(1,:) = [];
        minIntensities(1,:) = [];
        continue
    end
    %check values
    if(minIntensities(1)<minIntensities(currIdx))
        centroidToAppend = [currPointX,currPointY];
        centroidFromAppend = [centroidOrig(currIdx,1),centroidOrig(currIdx,2)];
    else
        centroidToAppend = [centroidOrig(currIdx,1),centroidOrig(currIdx,2)];
        centroidFromAppend = [currPointX,currPointY];
    end
    
    centroidTo = [centroidTo; centroidToAppend];
    centroidFrom = [centroidFrom; centroidFromAppend];
    
    centroidOrig(currIdx,:) = [];
    minIntensities(currIdx,:) = [];
    centroidOrig(1,:) = [];
    minIntensities(1,:) = [];
    
end

end

