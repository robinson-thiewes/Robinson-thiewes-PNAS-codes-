function af_new = flipGL(af, selLn)

af_new = af;
wl = length(selLn);




for i = 1:wl
    mn = size(af{selLn(i),7},2);
    
    for j =1:5
        af_new{selLn(i),j}(:,1) = mn - af{selLn(i),j}(:,1) + 1;
    end
    
    for j = [6 7 10 11]
        af_new{selLn(i),j} = flip(af{selLn(i),j},2);
    end
    
    for j =12
        af_new{selLn(i),j}(:,1) = mn - af{selLn(i),j}(:,1) + 1;
    end
    
    for j =14:15
        af_new{selLn(i),j}(:,1) = mn - af{selLn(i),j}(:,1) + 1;
    end
end

