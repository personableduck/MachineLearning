function [tailArray] = buildTail(currList, x_limit,y_limit, maxNodeCount, currTail, nodeLength, degreeSpan)
%Recursive function which builds tail projections for the tail specified in
%currList(currTail).  It returns an array of tails. 
%x_limit =  Max x index of image
%y_limit =  Max y index of image
%MaxLength =  Max number of nodes to grow a tail projection.
%currList = list or aray of current tail projections.
%currTail = the index of the current tail from current tail projections for
%which you want to build tail projections
%nodeLength = the length of the node to build.  (in pixels.  It was
%converted from microns in the script that calls this function.)
%degreeSpan = the angle span, and angle step size used to build tail projections.

currNodeCount = currList(currTail).nodeCount;
%Base condition check.  If the length of the projection has reached
%maximum length specified, stop building and return.
if(currNodeCount==maxNodeCount)
    tailArray = currList(currTail);
    return;
end


% currLength = currList(currTail).length;
% %Base condition check.  If the length of the projection has reached
% %maximum length specified, stop building and return.
% if(currLength== maxLength)
%     tailArray = currList(currTail);
%     return;
% end

orig_x = currList(currTail).x(currNodeCount);
orig_y = currList(currTail).y(currNodeCount);
orientation = currList(currTail).orientation;

if(orientation < 0)
    orientation = orientation+360;
end

degreeSpani = degreeSpan+orientation;
cPointList = zeros(size(length(degreeSpani),2));

%For each angle specified, build one node for the tial projection the size
%of node length specified.
for iter = 1:length(degreeSpani)
    
    m_orientation = degreeSpani(iter);
    rad = m_orientation*pi/180;
    x = orig_x + nodeLength*cos(rad);
    if (x<=0)
        x = 1;
    elseif(x>=x_limit)
        x = x_limit;
    end
    
    y = orig_y - nodeLength*sin(rad);
    if (y<=0)
        y = 1;
    elseif(y>=y_limit)
        y = y_limit;
    end
    
    cPointList(iter,1) = x;
    cPointList(iter,2) = y;
end

%take only unique points build
% cPointList = unique(cPointList,'rows');
numPoints = size(cPointList,1);

%copy current tail for as many points as discovered
[copyList] = replicateTail(currList(currTail),numPoints);
for i = 1:numPoints
	copyList(i,1).x(currNodeCount+1) = cPointList(i,1);
    copyList(i,1).y(currNodeCount+1) = cPointList(i,2);
    copyList(i,1).nodeCount = copyList(i,1).nodeCount+1;
    copyList(i,1).orientation = degreeSpani(i);
end

%for each tail made and orig tail passed, run function
tempArray = spermTail.empty(100,0,numPoints);
sizeArr = zeros(numPoints,1);
for i = 1:numPoints
    %grow by calling function again
    tailArray = buildTail(copyList,x_limit,y_limit,maxNodeCount,i, nodeLength, degreeSpan);
    tempArray(1:size(tailArray,1),1,i) = tailArray;
    sizeArr(i,1) = size(tailArray,1);
end

%After recursion, arrange all the tails into one array.
currIter = 1;
for i = 1:numPoints
    tailArray(currIter:currIter+sizeArr(i,1)-1,:) = tempArray(1:sizeArr(i,1),:,i);
    currIter = currIter+sizeArr(i,1);
end

end