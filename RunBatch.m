 %interval [ 0.3 0.4 0.5 ]; %%%% this is useful if you want unevenly spaced
%%%intervals
% interval1= [1 0.05 1.25 ] ; = 0 : 0.2 : 1  [0 : 0.2 : 1]
% interval2= [ 1 1 10  ]   ;
itvl_ex =  0.5;
itvl_int = 0.45;

siz1 = size(itvl_ex,2);
siz2 = size(itvl_int,2);


% % for both intron and exon channels
for i = 1:siz1
    for j = 1:siz2

        thresForExon =       itvl_ex(i)       ;
        thresForIntron =     itvl_int(j)    ;

        if i < 10
            pDir2 = strcat('batch results cond0',num2str(j+siz2*(i-1)),'_thrEx_', num2str(thresForExon), '_thrInt_', num2str(thresForIntron) );
        else
            pDir2 = strcat('batch results cond',num2str(j+siz2*(i-1)),'_thrEx_', num2str(thresForExon), '_thrInt_', num2str(thresForIntron) );
        end
        John_smFISH_v1_02_mpk1_v1_3(pDir2, thresForExon, thresForIntron);
    end
end


