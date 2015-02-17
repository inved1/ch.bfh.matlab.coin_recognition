threshold = 50;
array = [6712 7023 7510 7509 6718 7514 7509 6247];
sortedArray = sort(array);
nPerGroup = diff(find([1 (diff(sortedArray) > threshold) 1]));
groupArray = mat2cell(sortedArray,1,nPerGroup)