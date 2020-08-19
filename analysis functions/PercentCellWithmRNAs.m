function [] = PercentCellWithmRNAs(aff, rc, stLoc, iminfo)

%%% This function make two plots: (1) percentage of cell with mRNAs above
%%% basal level; and (2) number of mRNAs per cell, both along the germline
%%% from distal end. It gets four input variables: 
%%% 1st input variable: 'aff' or normalized smFISH analysis result variable 'af'.
%%% 2nd input variable: mRNA channel (e.g. rc = 1, being 1st RNA channel is
%%% for mRNA.
%%% 3rd input variable: sets the range to measure basal mRNA level (either
%%% background mRNA level or the threshold mRNA level that you want to set.
%%% The range is the relative location in micrometer from distal end. For
%%% example, stLoc = [ 40 60 ], means basal mRNA level will be measured
%%% within the region in the gonad from 40 to 60 micrometer from the distal end.
%%% 4th input variable: iminfo, information about the image
%%% Outputs are two plots as explained above.

%%% Usage: PercentCellWithmRNAs(aff, 1, [ 40 60 ], iminfo)

rc = rc*2-1;
 

crs=5;
xcr = 0:crs:60;
bins = length(xcr);
perATSpCR = zeros(size(aff,1), 13);
NmRNAperCell = zeros(size(aff,1), 13);
for i=1:size(aff,1)
    if stLoc(1) == stLoc(2)
        basalmRNAa = 0;
    else
    basalmRNAa = sum(aff{i,rc}(aff{i,rc}(:,1) >= stLoc(1)/iminfo(6) ...
        & aff{i,rc}(:,1) < stLoc(2)/iminfo(6), 11))/ (stLoc(2)-stLoc(1))*5;   
    end
    
    for j=1:bins     % 1 - 13th cell row
        xr = [crs*(j-1)/iminfo(6) crs*j/iminfo(6)];
        allnucinfo = aff{i,rc}(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2),10);
        allnucs = sum(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2));
        if allnucs == 0
            allnucs = 1;
        end
        perATSpCR(i,j) = sum(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2) ...
            & aff{i,rc}(:,10) > basalmRNAa) / allnucs * 100;
        NmRNAperCell(i,j) = mean(allnucinfo);
    end
end

meanATSpCR = mean(perATSpCR);
stdATSpCR = std(perATSpCR);
NmRNAperCellmean = mean(NmRNAperCell);
NmRNAperCellstd = std(NmRNAperCell);


%----------bargraph for % cell with mRNA above basal leve
figure

hold off
bar(xcr, meanATSpCR)
hold on
errorbar(xcr, meanATSpCR, stdATSpCR/sqrt(30), 'k.', 'linewidth', 2);
title('% cell with mRNA')
axis([ -4 64 0 100 ])
box on



%----------bargraph for # mRNA per cell
figure
hold off
bar(xcr, NmRNAperCellmean)
hold on
errorbar(xcr, NmRNAperCellmean, NmRNAperCellstd/sqrt(15), 'k.', 'linewidth', 2);
title('# of mRNA per cell')
axis([ -4 64 0 50 ])
box on

