function retval = checkExisting(numbers, nrNew, threshold)

    retval = false;
    
    min = nrNew - threshold;
    max = nrNew + threshold;

    for e = 1:size(numbers',1)
        
        nrSearched = numbers(e);
        
        if(nrSearched <= max)
            if(nrSearched >= min)
                retval = true;
            end
        end
    end
    
end