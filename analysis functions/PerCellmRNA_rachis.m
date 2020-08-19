function [NmRNAperCell, perMRNApCR,NmRNAperCellmean,NmRNAperCellstd] = PerCellmRNA(aff, mRNAch, y_max, stLoc, max_stLoc, iminfo)
%calculates the amount of mRNA per cell and in rough cell row groups for rachis
%stLoc and max_stLoc is range for basal mRNA (i.e. stLoc = 0 and max_stLoc = 10
%

rc = mRNAch;
rc = rc*2-1;
xc = mRNAch *2;

% crs = 4.55;       % cell row size (um)
% xum = 0:4.55:55;
% xcr = 0:1:60;


%%% METHOD 1: calculate basal mRNA level (95% CI of # mRNA per cell at proximal region 40-60um fr distal end)
% basalmRNA = zeros(size(aff,1),1);
% for i=1:size(aff,1)
%     tempCount = aff{i,1}(aff{i,1}(:,1) >= stLoc/iminfo(6) & aff{i,rc}(:,1) < max_stLoc/iminfo(6), 10);
%     sem = std(tempCount)/sqrt(length(tempCount));
%     ts = tinv([0.025 0.975], length(tempCount)-1);
%     CIup = mean(tempCount) + ts*sem;
%     basalmRNA(i) = CIup(2);    
% end


%%% METHOD 2: calculate top 5% cutoff.
basalmRNA = zeros(size(aff,1),1);
for i=1:size(aff,1)
    tempCount = aff{i,1}(aff{i,1}(:,1) >= stLoc/iminfo(6) & aff{i,rc}(:,1) < max_stLoc/iminfo(6), 9);
    tempC = sort(tempCount);
    if isempty(tempC)
        basalmRNA(i) = 0;
    else
        basalmRNA(i) = tempC(round(length(tempC)*97.5/100));   
    end
end

%  basalmRNAa = mean(basalmRNA);    


crs=5;
xcr = 0:crs:60;
bins = length(xcr);
perMRNApCR = zeros(size(aff,1), bins);
NmRNAperCell = zeros(size(aff,1), bins);
NmRNAperCell_aboveBasal = zeros(size(aff,1), bins);
mRNA_rachis = zeros(size(aff,1), bins);
for i=1:size(aff,1)
    %%% calculate basal mRNA level (mRNA # at proximal region 40-60um fr distal end)
    %%% L4: 45-65;    N2: 55-65;    JK5008: 40           lst-1: 55       lag-3/+:        .
%     stLoc = 40;
%     basalmRNAa = sum(aff{i,rc}(aff{i,rc}(:,1) >= stLoc/iminfo(6) ...
%         & aff{i,rc}(:,1) < 65/iminfo(6), 10))/ (65-stLoc)*5;    
    
    for j=1:bins     % 1 - 13th cell row
        xr = [crs*(j-1)/iminfo(6) crs*j/iminfo(6)];
        allnucinfo = aff{i,rc}(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2),10);
        allnucinfo2 = aff{i,rc}(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2) & basalmRNA(i) < 1 ,10);  %%% gonads with basal mRNA level < 20 are excluded from analysis. 
        allnucs = sum(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2));
        if allnucs == 0
            allnucs = 1;
        end
    end
        
%         rachis_mRNAs = aff{i,xc}(:,10) == 0; %make list of mRNAs in rachis based on Voronoi count with limit
%         cell_mRNAs = aff{i,xc}(:,10) == 1; %list of mRNA cell based on Voronoi count with limit
        mRNA_info = aff{i,xc}(:,:);
%         rachis = mRNA_info;
        cellmRNA = mRNA_info;
%         rachis(cell_mRNAs,:) = []; %removes mRNA in cell and leaves only rachis mRNAs
%         cellmRNA(rachis_mRNAs,:) = []; %removes rachis mRNA and leaves cell mRNAs
        
        
        perMRNApCR(i,j) = sum(cellmRNA(:,1) >= xr(1) & cellmRNA(:,1) < xr(2) ...
            & cellmRNA(:,10) > basalmRNA(i)) / allnucs * 100;
        NmRNAperCell(i,j) = mean(allnucinfo);
        
        NmRNAperCell_aboveBasal(i,j) = mean(allnucinfo2(allnucinfo2 > basalmRNA(i)));
        
%         mRNA_rachis(i,j) = sum(rachis(:,1) >= xr(1) &rachis(:,1) < xr(2)) ;
        
%     end
end

perMRNApCR(isnan(perMRNApCR)) = 0;
NmRNAperCell(isnan(NmRNAperCell)) = 0;
NmRNAperCell_aboveBasal(isnan(NmRNAperCell_aboveBasal)) = 0;

% mean_mRNA_rachis = mean(rachis);
% std_mRNA_rachis = std(rachis);

meanMRNApCR = mean(perMRNApCR);
stdMRNApCR = std(perMRNApCR);
NmRNAperCellmean = mean(NmRNAperCell);
NmRNAperCellstd = std(NmRNAperCell);

NmRNAperCell_aBasalmean = zeros(1, size(NmRNAperCell_aboveBasal,2));
NmRNAperCell_aBasalstd = zeros(1, size(NmRNAperCell_aboveBasal,2));
for i = 1:size(NmRNAperCell_aboveBasal,2)
    NmRNAperCell_aBasalmean(i) = mean(NmRNAperCell_aboveBasal(NmRNAperCell_aboveBasal(:,i)>0,i));
    NmRNAperCell_aBasalstd(i) = std(NmRNAperCell_aboveBasal(NmRNAperCell_aboveBasal(:,i)>0,i));
end



%----------bargraph for % cell with mRNA above basal leve
figure('pos', [300 200 350 500])

hold off
h1=bar(xcr+crs/2, meanMRNApCR);
h1.FaceColor = [1 0 .5];
hold on
errorbar(xcr+crs/2, meanMRNApCR, stdMRNApCR/sqrt(size(perMRNApCR,1)), 'k.', 'linewidth', 2);
axis([ -2 64 0 100 ])
xticks(0:5:60)
box on
ylabel('Percent os cells with \itmpk-1\rm mRNA' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);

% %----------bargraph for #mRNA in rachis 
% figure('pos', [300 200 350 500])
% 
% hold off
% h1=bar(xcr+crs/2, mean_mRNA_rachis);
% h1.FaceColor = [1 0 .5];
% hold on
% errorbar(xcr+crs/2, mean_mRNA_rachis, std_mRNA_rachis/sqrt(size(rachis,1)), 'k.', 'linewidth', 2);
% axis([ -2 64 0 100 ])
% xticks(0:5:60)
% box on
% ylabel('Percent os cells with \itmpk-1\rm mRNA' , 'fontsize',15);
% xlabel('\itum\rm from distal end', 'fontsize',15);


%----------bargraph for # mRNA per cell
figure('pos', [300 200 350 500])
hold off
h2=bar(xcr+crs/2, NmRNAperCellmean);
h2.FaceColor = [1 0 .5];
hold on
xticks(0:5:60)
errorbar(xcr+crs/2, NmRNAperCellmean, NmRNAperCellstd/sqrt(size(NmRNAperCell,1)), 'k.', 'linewidth', 2);
axis([ -0 60 0 y_max ])
ylabel('Number \itmpk-1\rm mRNA per cell' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);

%---------boxplot for # mRNA per cell
% figure('pos', [300 200 350 500])
% boxplot(NmRNAperCell, 'Notch', 'on', 'Whisker', 1)
% hold on
% ylabel('Number \itmpk-1\rm mRNA per cell' , 'fontsize',15);
% xlabel('cell rows from distal end', 'fontsize',15);




% %-------(optional) bargraph for # mRNA per cell for only cell with mRNA # above basal level.
figure('pos', [300 200 350 500])
hold off
h3=bar(xcr+crs/2, NmRNAperCell_aBasalmean);
h3.FaceColor = [1 0 .5];
hold on
xticks(0:10:100)
errorbar(xcr+crs/2, NmRNAperCell_aBasalmean, NmRNAperCell_aBasalstd/sqrt(size(NmRNAperCell_aboveBasal,1)), 'k.', 'linewidth', 2);
axis([ -4 64 0 y_max ])
ylabel({'(3) # \itmpk-1\rm mRNA per cell';'(only cells with mRNA above basal level)'} , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);





