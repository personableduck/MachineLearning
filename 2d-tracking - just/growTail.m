initList = spermTail.empty(1,0);
orientation_found = tailArray(maxIdx,1).orientation;


initList(1,1).x = tailArray(maxIdx,1).x(end-1);
initList(1,1).y = tailArray(maxIdx,1).y(end-1);
initList(1,1).x(2) = tailArray(maxIdx,1).x(end);
initList(1,1).y(2) = tailArray(maxIdx,1).y(end);
initList(1,1).nodeCount = overlap;
% initList(1,1).length = 2;
initList(1,1).orientation = tailArray(maxIdx,1).orientation;
currTail = 1;

tailArray = buildTail(initList, x_limit, y_limit, maxNodeCount, currTail, nodeLength, degreeSpan);
