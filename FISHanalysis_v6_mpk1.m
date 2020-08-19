%% FISHAnalysis_v6_mpk-1 %% 
%Use to analyse smFISH images for mpk-1
%Name code:
%N2HR1= N2 hermaphrodite replicate 1
%N2HR2= N2 hermaphrodite replicate 2
%N2MR1= N2 male replicate 1
%N2MR2= N2 male replicate 2
%DTHR1= q1083 (dual tagged) hermphrodite replicate 1
%DTHR2= q1083 (dual tagged) hermaphrodite replicate 2
%DTMR1= q1083 (dual tagged) male replicate 1
%DTMR2= q1083 (dual tagged) male replicate 2
%ICHR1= intron control hermaphrodite replicate 1
%ICHR2= intron control hermpahridote repilcate 2
 
%% Flip germlines, must be done manually%% 
%axis orientation: DTC on left
%only need to flip those that are missed by the first flip (incorrect direction when checking RNA detection)
%flipGL(name, [#s to be flipped])

af_test = flipGL(af, [4 5]);
% if exist('af_N2HR1','var') 
%     af_N2HR1 = flipGL(af_N2HR1, [ 26 27 37 53]);
% end
% 
% if exist('af_N2HR2','var') 
%     af_N2HR2 = flipGL(af_N2HR2, [ 26 27 37 53]);
% end
% 
% if exist('af_N2MR1','var') 
%     af_N2MR1 = flipGL(af_N2MR1, [ 26 27 37 53]);
% end
% 
% if exist('af_N2MR2','var') 
%     af_N2MR2 = flipGL(af_N2MR2, [ 26 27 37 53]);
% end
% 
% if exist('af_DTHR1','var') 
%     af_DTHR1 = flipGL(af_DTHR1, [ 26 27 37 53]);
% end
% 
% if exist('af_DTHR2','var') 
%     af_DTHR2 = flipGL(af_DTHR2, [ 26 27 37 53]);
% end
% 
% if exist('af_DTMR1','var') 
%     af_DTMR1 = flipGL(af_DTMR1, [ 26 27 37 53]);
% end
% 
% if exist('af_DTMR2','var') 
%     af_DTMR2 = flipGL(af_DTMR2, [ 26 27 37 53]);
% end

%% rotate germline
% 
% [af_new1] = rotate(af,1,210);
% [af_new2] = rotate(af,2,200);
% [af_new4] = rotate(af,4,220);
%  [af_new6] = rotate(af,6,245);
%% Process image data before further analyses: Normalization %%
%%% ============= Align coordinates in 'aff' =======================

%----- METHOD 1: align by left end of the first nucleus (0 on x-axis) -----

% For 2nd input variable, designate the RNA channels.
% Put   '1' for intron,     '2' for exon,      '3' for else.
% For 3rd input, list n-th images to exclude from normalization.
[aff, mrnaNucAll, mrnaAll, atsNucAll, atsAll] = normaf(af, [1 2], []); %, intron_values_all, mean_intron_values);

[aff, mrnaNucAll, mrnaAll, ~, ~] = normaf(af, [1 2], [7]) %, intron_values_all, mean_intron_values);

%For N2 hermaphrodite replicate 1 
if exist('af_N2HR1','var')&& exist('iminfo_N2HR1','var')
    [aff_N2HR1, mrnaNucAll_N2HR1, mrnaAll_N2HR1, atsNucAll_N2HR1, atsAll_N2HR1] = normaf(af_N2HR1, [ 1  2 ], []); %af, [channel order], [numbers to skip]
end 

%For N2 hermaphrodite replicate 2
if exist('af_N2HR2','var') && exist('iminfo_N2HR2','var')
    [aff_N2HR2, mrnaNucAll_N2HR2, mrnaAll_N2HR2, atsNucAll_N2HR2, atsAll_N2HR2] = normaf(af_N2HR2, [ 1  2 ], [3]);
end

%For N2 male replicate 1
if exist('af_N2MR1', 'var') && exist('iminfo_N2MR1','var')
    [aff_N2MR1, mrnaNucAll_N2MR1, mrnaAll_N2MR1, atsNucAll_N2MR1, atsAll_N2MR1] = normaf(af_N2MR1, [ 1  2 ], []);
end

%For N2 male replicate 2
if exist('af_N2MR2','var') && exist('iminfo_N2MR2','var')
    [aff_N2MR2, mrnaNucAll_N2MR2, mrnaAll_N2MR2, atsNucAll_N2MR2, atsAll_N2MR2] = normaf(af_N2MR2, [ 1  2 ], []);
end

%For dual tagged hermaphrodite repicate 1
if exist('af_DTHR1','var') && exist('iminfo_DTHR1','var')
    [aff_DTHR1, mrnaNucAll_DTHR1, mrnaAll_DTHR1, atsNucAll_DTHR1, atsAll_DTHR1] = normaf(af_DTHR1, [ 1  2 ], []);
end

%For dual tagged hermpahrodite replicate 2
if exist('af_DTHR1','var') && exist('iminfo_DTHR2','var')
    [aff_DTHR2, mrnaNucAll_DTHR2, mrnaAll_DTHR2, atsNucAll_DTHR2, atsAll_DTHR2] = normaf(af_DTHR2, [ 1  2 ], []);
end

%For dual tagged male replicate 1
if exist('af_DTMR1','var') && exist('iminfo_DTMR1','var')
    [aff_DTMR1, mrnaNucAll_DTMR1, mrnaAll_DTMR1, atsNucAll_DTMR1, atsAll_DTMR1] = normaf(af_DTMR1, [ 1  2 ], []);
end

%For dual tagged male replicate 2
if exist('afDTMR2','var') && exist('iminfo_DTMR2','var')
    [aff_DTMR2, mrnaNucAll_DTMR2, mrnaAll_DTMR2, atsNucAll_DTMR2, atsAll_DTMR2] = normaf(af_DTMR2, [ 1  2 ], []);
end

%% Number of cells per cell row
%[perRowNuclei, mean_perRowNuclei, std_perRowNuclei]=NumCells(aff, ATSch, iminfo)

[perRowNuclei, mean_perRowNuclei, std_perRowNuclei]=NumCells(aff, 1, iminfo)

%% Intenisty of a single mRNA

for i=1:size(aff,1)
    mrnas = aff{i,2}(:,:);
    cyto_mRNA = mrnaAll;
    nuc_calledmRNA = mrnaAll(:,12) == 0;
    cyto_mRNA(nuc_calledmRNA,:) = [];
    cyto_mRNA_all{i} = cyto_mRNA;
    mean_cyto_mRNA{i} = mean(cyto_mRNA(:,6));
    std_cyto_mRNA{i} = std(cyto_mRNA(:,6));
end

means_cyto = cell2mat(mean_cyto_mRNA);
mrnaAllp = mrnaAll;
cyto_mRNA = mrnaAllp;
nuc_mRNA = mrnaAllp;
nuc_calledmRNA = mrnaAll(:,12) == 0;
cyto_calledmRNA = mrnaAll(:,12) > 0;
cyto_mRNA(nuc_calledmRNA,:) = [];
nuc_mRNA(cyto_calledmRNA,:) = [];
cyto_mRNA(cyto_mRNA(:,6) > mean(cyto_mRNA(:,6))*4,5) = cyto_mRNA(cyto_mRNA(:,6) > mean(cyto_mRNA(:,6))*4,6)/mean(cyto_mRNA(:,6))*2;



figure
hist(cyto_mRNA(:,6)./cyto_mRNA(:,5), 30)
title('Signal intensity of mRNA (\ita.u.\rm)', 'fontsize',15);
axis([ 0 60 0 5000])
CoefOfVariation = std(cyto_mRNA(:,6)./cyto_mRNA(:,5))/mean(cyto_mRNA(:,6)./cyto_mRNA(:,5))
% 
% figure
% hist(nuc_mRNA(:,6)./nuc_mRNA(:,5), 30)
% title('Signal intensity of mRNA (\ita.u.\rm)', 'fontsize',15);
% axis([ 0 60 0 1*1e4 ])
% CoefOfVariation = std(nuc_mRNA(:,6)./nuc_mRNA(:,5))/mean(nuc_mRNA(:,6)./nuc_mRNA(:,5))

%% exon ATS
for i=1:size(aff,1) %i = # of images
    [~,ia,ib] = intersect(aff{i,2}(:,8), aff{i,4}(:,9), 'rows', 'stable'); %find intersection between nuc ID listed in intron ATS info and exon nuc info
    exon_ATS = aff{i,2}(ia,:); %index interesected RNA IDs
    exons_ATS{i} = exon_ATS;
    mean_exons_ATS{i} = mean(exon_ATS(:,6));
    std_exons_ATS{i} = std(exon_ATS(:,6));
end
    means_ATS_exon = cell2mat(mean_exons_ATS);
    test_exons = cell2mat(exons_ATS');
    test_exons(test_exons(:,6) > mean(test_exons(:,6))*4,5) = test_exons(test_exons(:,6) > mean(test_exons(:,6))*4,6)/mean(test_exons(:,6))*2;

    figure
    hist(test_exons(:,6)./test_exons(:,5), 30)
    title('Signal intensity of exon ATS (\ita.u.\rm)', 'fontsize',15);
    axis([ 0 60 0 1000])
    CoefOfVariation = std(test_exons(:,6)./test_exons(:,5))/mean(test_exons(:,6)./test_exons(:,5))

    [~,ia,ib] = intersect(aff{1,4}(:,8), aff{1,2}(:,9), 'rows', 'stable');
    test_exon = aff{1,4}(ia,:);
    
    [~,ia,ib] = intersect(aff{1,4}(:,9), aff{1,2}(:,8), 'rows', 'stable');
    test_exon = aff{1,4}(ib,:);
        
%% Intensity of ATS (intron)
%intron
atsAllp = atsAll;
atsAllp(atsAllp(:,6) > mean(atsAllp(:,6))*4,5) = atsAllp(atsAllp(:,6) > mean(atsAllp(:,6))*4,6)/mean(atsAllp(:,6))*2;

figure
hist(atsAllp(:,6)./atsAllp(:,5), 30)
title('Signal intensity of intron ATS (\ita.u.\rm)', 'fontsize',15);
axis([ 0 60 0 2000])
CoefOfVariation = std(atsAll(:,6)./atsAll(:,5))/mean(atsAll(:,6)./atsAll(:,5))

%----------- % of mRNA spots were detected as 1 mRNA.
x= sum(mrnaAll(mrnaAll(:,5) == 1,4))/size(mrnaAll,1)*100
%% ==================== % cell with ATS in space ==============================
%%% This code draws plot for % of cells with ATS against distance (um) from distal end.

[perATSpCR, meanATSpCR, stdATSpCR]=PerCellATS(aff, 1, iminfo);

%combined hermaphrodite data
[perATSpCR_N2HRC, meanATSpCR_N2HRC, stdATSpCR_N2HRC]=PerCellATS(N2HRC, 1, iminfo);

%combined male data
[perATSpCR_N2MRC, meanATSpCR_N2MRC, stdATSpCR_N2MRC]=PerCellATS(N2MRC, 1, iminfo);

%graph of intron and exon detected ATS

figure('pos', [300 200 350 500])
hold off
plot(xcr+crs/2, mpk_meanATSpCR, xcr+crs/2, lag3_meanATSpCR, xcr+crs/2, lag1_meanATSpCR, xcr+crs/2, sygl_meanATSpCR);
axis([0 60 0 100])
plot(xcr+crs/2, mpk_mean_mRNA_production, xcr+crs/2, lag3_mean_mRNA_production, xcr+crs/2, lag1_mean_mRNA_production, xcr+crs/2, sygl_mean_mRNA_production);
axis([0 60 0 15])
hold on
pause(0.1)
xticks(0:5:60);
axis([ -2 60 0 100 ]);
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
legend('Intron detected ATS', 'Exon detected ATS', 'location', 'northeast')
ylabel('% cell with \itmpk-1b\rm ATS' , 'fontsize',15);
xlabel('\itum\rm in window', 'fontsize',15);

%% hermaprodite vs male stats and graphs

crs=5;
xcr = 0:crs:60;
[h_N2HM_ATS, p_N2HM_ATS] = ttest2(perATSpCR_N2HRC,perATSpCR_N2MRC);
figure
hold on
hbar=bar(xcr,[meanATSpCR_N2HRC',meanATSpCR_N2MRC']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,meanATSpCR_N2HRC',stdATSpCR_N2HRC'/sqrt(size(perATSpCR_N2HRC,1)), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,meanATSpCR_N2MRC',stdATSpCR_N2MRC'/sqrt(size(perATSpCR_N2MRC,1)), 'k.', 'linewidth', 2);
xticks(0:5:60);
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
ylabel('%cells with \itmpk-1\rm NTS', 'fontsize', 15);
xlabel('\itum\rm in window', 'fontsize', 15);
legend('Wild type hermaphrodite', 'Wild type male', 'location', 'northeast');
hold off

%% 2-tailed unpaired T test and Wilcoxon test comparing replicates for cell with ATS in space
%use to compare if replicates are significally different from each other for ATS in space
%[h,p]=ttest2(x,y)   x = replicate 1 and y = replicate 2
%[p,h]=ranksum(x,y) Wilcoxon test  
% h = 0, fail to reject the null hypothesis, h = 1 reject null
%N2HRC= concantenated replicates 1 and 2 hermaphrodite

%For N2 hermpahrodite replicates
if exist('perATSpCR_N2HR1', 'var') && exist('perATSpCR_N2HR2', 'var')
    [h_N2HRATS, p_N2HRATS]= ttest2(perATSpCR_N2HR1, perATSpCR_N2HR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
    if h_N2HRATS(1)==0 %combine dataset only if replicates not significantly different
        N2HRC_ATSpercent = vertcat(perATSpCR_N2HR1, perATSpCR_N2HR2);
    end
end

%For N2 male replicates
if exist('perATSpCR_N2MR1', 'var') && exist('perATSpCR_N2MR2', 'var')
    [h_N2MRATS, p_N2MRATS]= ttest2(perATSpCR_N2MR1, perATSpCR_N2MR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
    if h_N2MRATS==0 %combine dataset only if replicates not significantly different
        N2MRC_ATSpercent = vertcat(perATSpCR_N2MR1, perATSpCR_N2MR2);
    end 
end

%For dual tagged hermaphrodite replicates
if exist('perATSpCR_DTHR1', 'var') && exist('perATSpCR_DTHR2', 'var')
    [h_DTHATS, p_DTHATS]= ttest2(perATSpCR_DTHR1, perATSpCR_DTHR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
    %if h_DTHATS==0 %combine dataset only if replicates not significantly different
    DTHRC_ATSpercent = vertcat(perATSpCR_DTHR1, perATSpCR_DTHR2);
    %end 
end

%For dual tagged male replicates
if exist('perATSpCR_DTHM1', 'var') && exist('perATSpCR_DTHM2', 'var')
    [h_DTMATS, p_DTMATS]= ttest2(perATSpCR_DTHM1, perATSpCR_DTHRM2);
    h_DTMATS(13,:)=[]; %put in line to delete 13th column because no values
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
    if h_DTMRATS==0 %combine dataset only if replicates not significantly different
        DTMRC_ATSpercent = concat(perATSpCR_DTMR1, perATSpCR_DTMR2);
    end 
end



%% Comparison between genotypes for ATS in space
% 2 tailed unpaired t-test
%put in lines to make figure

%N2 vs dual tagged hermahrodites
if exist('N2HRC_ATSpercent', 'var') && exist('DTHRC_ATSpercent', 'var') 
    [h_N2vsDTH, p_N2vsDTH]= ttest2(N2HRC_ATSpercent, DTHRC_ATSpercent);
    meanN2HRC_ATSpercent = mean(N2HRC_ATSpercent);
    meanDTHRC_ATSpercent = mean(DTHRC_ATSpercent);
    stdN2HRC_ATSpercent = std(N2HRC_ATSpercent);
    stdDTHRC_ATSpercent = std(DTHRC_ATSpercent);
    [h_N2vsDTH_ATSpercent, p_N2vsDTH_ATSpercent]=ttest2(meanN2HRC_ATSpercent,meanDTHRC_ATSpercent);
    figure
    hold on
    hbar=bar(xcr,[meanN2HRC_ATSpercent']);
    pause(0.1)
    errorbar(xcr+hbar(1).XOffset,meanN2HRC_ATSpercent',stdN2HRC_ATSpercent'/sqrt(50), 'k.', 'linewidth', 2);
%     errorbar(xcr+hbar(2).XOffset,meanDTHRC_ATSpercent',stdDTHRC_ATSpercent'/sqrt(50), 'k.', 'linewidth', 2);
    legend('N2 hermaprodite', 'location', 'northeast');
    xticks(0:10:100)
    hold off
end
figure
hold on
hbar = bar(xcr,meanN2HRC_ATSpercent)
pause(0.1)
errorbar(xcr+hbar(1).XOffset, meanN2HRC_ATSpercent, stdN2HRC_ATSpercent/sqrt(50), 'k.', 'linewidth', 2)
legend('N2 hermaphrodite', 'location', 'northeast');
xticks(0:10:100)
hold off



%For N2 vs dual tagged males 
if exist('N2MRC_ATSpercent', 'var') && exist('DTMRC_ATSpercent', 'var') 
    [h_N2vsDTM, p_N2vsDTM]= ttest2(N2MRC_ATSpercent, DTMRC_ATSpercent);
end

%% Comparison between sexes for ATS in space
%2 tailed unpaired t-test

%N2 hermaphrodite vs N2 males
if exist('N2HRC_ATSpercent', 'var') && exist('N2MRC_ATSpercent', 'var') 
    [h_N2HvM, p_N2HvM]= ttest2(N2HRC_ATSpercent, N2MRC_ATSpercent);
end

%For dual tagged hermaphrodites vs dual tagged males
if exist('DTHRC_ATSpercent', 'var') && exist('DTMRC_ATSpercent', 'var') 
    [h_DTHvM, p_DTHvM]= ttest2(DTHRC_ATSpercent, DTMRC_ATSpercent);
end


    
    


%% =================== Trx activity (summed ATS intensity) in the germline =================
%Plots radisus of nucleus and DAPI intensity
%how to get correlation? trend line how good of fit?
%lots of zeros at bottom

ATSintSpace(aff, 2, iminfo, 2000, atsNucAll);


%% ================= % and # cell with mRNA (above basal level) ==================
%input (aff, mRNA channel, y_max, stLoc, max_stLoc, iminfo)

%at moment, basal is any cell with more than 1 mRNA

%test code
[NmRNAperCell, perMRNApCR,NmRNAperCellmean,NmRNAperCellstd,TotalmRNAperCell]=PerCellmRNA(aff, 2, 15, 0, 1, iminfo);


%N2 hermaphrodite
[NmRNAperCell_N2HRC, perMRNApCR_N2HRC,NmRNAperCellmean_N2HRC,NmRNAperCellstd_N2HRC]=PerCellmRNA(N2HRC, 2, 20, 0, 1, iminfo);

%N2 male
[NmRNAperCell_N2MRC, perMRNApCR_N2MRC,NmRNAperCellmean_N2MRC,NmRNAperCellstd_N2MRC]=PerCellmRNA(N2MRC, 2, 20, 0, 1, iminfo);

%% Hermaphrodite vs male, stats and graphs
[h_NmRNAperCell_combined, p_NmRNAperCell_combined] = ttest2(NmRNAperCell_N2HRC,NmRNAperCell_N2MRC);

figure
hold on
hbar=bar(xcr,[NmRNAperCellmean_N2HRC',NmRNAperCellmean_N2MRC']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,NmRNAperCellmean_N2HRC',NmRNAperCellstd_N2HRC'/sqrt(size(NmRNAperCell_N2HRC,1)), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,NmRNAperCellmean_N2MRC',NmRNAperCellstd_N2MRC'/sqrt(size(NmRNAperCell_N2MRC,1)), 'k.', 'linewidth', 2);
legend('N2 hermaprodite', 'location', 'northeast');
xticks(0:5:60);
axis([-4 64 0 25])
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
ylabel('Number \itmpk-1\rm RNA per cell', 'fontsize', 15);
xlabel('\itum\rm in window', 'fontsize', 15);
legend('Wild type hermaphrodite', 'Wild type male', 'location', 'northeast');
hold off
%% t-test to compare between replicates

%use to compare if replicates are significally different from each other for # mRNA per cell
%[h,p]=ttest2(x,y)   x = replicate 1 and y = replicate 2
%[p,h]=ranksum(x,y) Wilcoxon test  
% h = 0, fail to reject the null hypothesis, h = 1 reject null
%N2HRC= concantenated replicates 1 and 2 hermaphrodite

%For N2 hermpahrodite replicates
if exist('NmRNAperCell_N2HR1', 'var') && exist('NmRNAperCell_N2HR2', 'var')
    [h_N2HRNmRNA, p_N2HRNmRNA]= ttest2(NmRNAperCell_N2HR1, NmRNAperCell_N2HR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
    N2HRC_NmRNA = vertcat(NmRNAperCell_N2HR1, NmRNAperCell_N2HR2);
end

%For N2 male replicates
if exist('NmRNAperCell_N2MR1', 'var') && exist('NmRNAperCell_N2HR2', 'var')
    [h_N2MRNmRNA, p_N2MRNmRNA]= ttest2(NmRNAperCell_N2MR1, NmRNAperCell_N2MR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
    if h_N2MRNmRNA==0 %combine dataset only if replicates not significantly different
        N2MRC_NmRNA = concat(NmRNAperCell_N2MR1, NmRNAperCell_N2MR2);
    end 
end

%For Dual Tagged hermpahrodite replicates
if exist('NmRNAperCell_DTHR1', 'var') && exist('NmRNAperCell_DTHR2', 'var')
    [h_DTHRNmRNA, p_DTHRNmRNA]= ttest2(NmRNAperCell_DTHR1, NmRNAperCell_DTHR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
     DTHRC_NmRNA = vertcat(NmRNAperCell_DTHR1, NmRNAperCell_DTHR2);
end

%For N2 hermpahrodite replicates
if exist('NmRNAperCell_N2HR1', 'var') && exist('NmRNAperCell_N2HR2', 'var')
    [h_N2HRNmRNA, p_N2HRNmRNA]= ttest2(NmRNAperCell_N2HR1, NmRNAperCell_N2HR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)
    if h_N2HRNmRNA==0 %combine dataset only if replicates not significantly different
        N2HRC_NmRNA = concat(NmRNAperCell_N2HR1, NmRNAperCell_N2HR2);
    end 
end

%% t-Test between genotypes

%N2 hermaphrodites vs dual tagged hermaphrodites
if exist('N2HRC_NmRNA','var') && exist('DTHRC_NmRNA', 'var')
    [h_N2vsDT_NmRNA, p_N2vsDT_NmRNA]=ttest2(N2HRC_NmRNA, DTHRC_NmRNA);
    mean_N2HRC_NmRNA = mean(N2HRC_NmRNA);
    std_N2HRC_NmRNA= std(N2HRC_NmRNA);
    mean_DTHRC_NmRNA=mean(DTHRC_NmRNA);
    std_DTHRC_NmRNA=std(DTHRC_NmRNA);
    figure
    hold on
    hbar=bar(xcr,[mean_N2HRC_NmRNA']);
    pause(0.1)
    errorbar(xcr+hbar(1).XOffset,mean_N2HRC_NmRNA',std_N2HRC_NmRNA'/sqrt(50), 'k.', 'linewidth', 2);
%     errorbar(xcr+hbar(2).XOffset,mean_DTHRC_NmRNA',std_DTHRC_NmRNA'/sqrt(25), 'k.', 'linewidth', 2);
    legend('N2 hermaprodite', 'location', 'northwest');
    xticks(0:10:100)
    hold off
end

figure
hold on
hbar=bar(xcr, mean_N2HRC_NmRNA')
pause(0.1)
errorbar(xcr+hbar(1).XOffset,mean_N2HRC_NmRNA',std_N2HRC_NmRNA'/sqrt(50), 'k.', 'linewidth', 2);
legend('N2 hermaphrodite', 'location', 'northwest');
xticks(0:10:100)
hold off

%% Production of mRNA per ATS over distance
%production(aff, iminfo,Exon, y)

[mRNA_production, mean_mRNA_production, std_mRNA_production]=production(aff, iminfo, 2, 15);

%N2 hermaphrodite
[mRNA_production_N2HRC, mean_mRNA_production_N2HRC, std_mRNA_production_N2HRC]=production(N2HRC, iminfo, 2, 15);

%N2 male
[mRNA_production_N2MRC, mean_mRNA_production_N2MRC, std_mRNA_production_N2MRC]=production(N2MRC, iminfo, 2, 15);
%% Hermaphrodite vs male, stats and graphs
[h_production_combined, p_production_combined] = ttest2(mRNA_production_N2HRC, mRNA_production_N2MRC);

figure
hold on
hbar=bar(xcr,[mean_mRNA_production_N2HRC',mean_mRNA_production_N2MRC']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,mean_mRNA_production_N2HRC',std_mRNA_production_N2HRC'/sqrt(size(mRNA_production_N2HRC,1)), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,mean_mRNA_production_N2MRC',std_mRNA_production_N2MRC'/sqrt(size(mRNA_production_N2MRC,1)), 'k.', 'linewidth', 2);
legend('N2 hermaprodite', 'location', 'northeast');
xticks(0:5:60);
axis([-4 64 0 15])
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
ylabel('Number \itmpk-1\rm RNA produced per ATS', 'fontsize', 15);
xlabel('\itum\rm in window', 'fontsize', 15);
legend('Wild type hermaphrodite', 'Wild type male', 'location', 'northeast');
hold off
%% Number of mRNA in rachis
% rachis(aff, Exon, iminfo, y)

[rachis_mRNA, rachis_mRNA_mean, rachis_mRNA_std]=Rachis(aff, 2, iminfo, 350);

%N2 hermaphrodite
[rachis_mRNA_N2HRC, rachis_mRNA_mean_N2HRC, rachis_mRNA_std_N2HRC]=Rachis(N2HRC, 2, iminfo, 50);

%N2 male
[rachis_mRNA_N2MRC, rachis_mRNA_mean_N2MRC, rachis_mRNA_std_N2MRC]=Rachis(N2MRC, 2, iminfo, 50);

%% hermaphrodite vs males, stats and graphs
[h_rachis_combined, p_rachis_combined] = ttest2(rachis_mRNA_N2HRC,rachis_mRNA_N2MRC);

figure
hold on
hbar=bar(xcr,[rachis_mRNA_mean_N2HRC',rachis_mRNA_mean_N2MRC']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,rachis_mRNA_mean_N2HRC',rachis_mRNA_std_N2HRC'/sqrt(size(rachis_mRNA_N2HRC,1)), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,rachis_mRNA_mean_N2MRC',rachis_mRNA_std_N2MRC'/sqrt(size(rachis_mRNA_N2MRC,1)), 'k.', 'linewidth', 2);
legend('N2 hermaprodite', 'location', 'northeast');
xticks(0:5:60);
axis([-4 64 0 250])
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
ylabel('Number \itmpk-1\rm RNA in the rachis', 'fontsize', 15);
xlabel('\itum\rm in window', 'fontsize', 15);
legend('Wild type hermaphrodite', 'Wild type male', 'location', 'northeast');
hold off

%% Correlation between ATS intron and exon intensity

overlapIntensity_correlation(atsAll_q1040, iminfo);

%% Burst analysis
[overlap_intronIntensity, overlap_exonIntensity, iATS_intensity]=Burst_intensity(atsAll);

group = [ones(size(A)); 2*ones(size(B)); 3*ones(size(C)); 4*ones(size(D))];
figure
boxplot([A; B; C; D], group)
ylim([0 50])
set(gca, 'XTickLabel', {'wildtype', '5 prime intron deletion', 'RNAP II peak deletion', '3 prime intron deletion'})
title('iATS intensity')


%% %intron site detection per ATS
[oneATS_count, twoATS_count, threeATS_count, fourATS_count] = atsCounts(aff, iminfo, 1);

[oneATS_count_male, twoATS_count_male, threeATS_count_male, fourATS_count_male] = atsCounts(N2MRC, iminfo, 1);


%% Percent overlap for ATS over distance and correlation tests
%function Overlap(aff, ATSCh, iminfo)
%for N2 hermaprodite replicate 1
%
[perATSoverlap,meanATSoverlap,stdATSoverlap,perATSIntron,meanATSIntron,stdATSIntron]=Overlap_v2(aff, 1, iminfo);

%N2 hermaphrodite
[perATSoverlap_N2HRC,meanATSoverlap_N2HRC,stdATSoverlap_N2HRC,perATSIntron_N2HRC,meanATSIntron_N2HRC,stdATSIntron_N2HRC]=Overlap_v2(N2HRC, 1, iminfo);

%N2 male
[perATSoverlap_N2MRC,meanATSoverlap_N2MRC,stdATSoverlap_N2MRC,perATSIntron_N2MRC,meanATSIntron_N2MRC,stdATSIntron_N2MRC]=Overlap_v2(N2MRC, 1, iminfo);

%% Hermaphrodite vs male, stats and graphs
[h_intron_combined, p_intron_combined] = ttest2(perATSIntron_N2HRC,perATSIntron_N2MRC);
[h_overlap_combined, p_overlap_combined] = ttest2(perATSoverlap_N2HRC,perATSoverlap_N2MRC);

%intron only figure
figure
hold on
hbar=bar(xcr,[meanATSIntron_N2HRC',meanATSIntron_N2MRC']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,meanATSIntron_N2HRC',stdATSIntron_N2HRC'/sqrt(size(perATSIntron_N2HRC,1)), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,meanATSIntron_N2MRC',stdATSIntron_N2MRC'/sqrt(size(perATSIntron_N2MRC,1)), 'k.', 'linewidth', 2);
xticks(0:5:60);
axis([-4 64 0 100])
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
ylabel('% intron only signal for \itmpk-1\rm ATS', 'fontsize', 15);
xlabel('\itum\rm in window', 'fontsize', 15);
legend('Wild type hermaphrodite', 'Wild type male', 'location', 'northeast');
hold off

%overlap figure
figure
hold on
hbar=bar(xcr,[meanATSoverlap_N2HRC',meanATSoverlap_N2MRC']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,meanATSoverlap_N2HRC',stdATSoverlap_N2HRC'/sqrt(size(perATSoverlap_N2HRC,1)), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,meanATSoverlap_N2MRC',stdATSoverlap_N2MRC'/sqrt(size(perATSoverlap_N2MRC,1)), 'k.', 'linewidth', 2);
xticks(0:5:60);
axis([-4 64 0 100])
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
ylabel('% intron and exon signal overlap at \itmpk-1\rm ATS', 'fontsize', 15);
xlabel('\itum\rm from distal end', 'fontsize', 15);
legend('Wild type hermaphrodite', 'Wild type male', 'location', 'northeast');
hold off

%% 2-tailed unpaired T test and Wilcoxon test comparing replicates for cell with ATS in space
%use to compare if replicates are significally different from each other for ATS in space
%[h,p]=ttest2(x,y)   x = replicate 1 and y = replicate 2
%[p,h]=ranksum(x,y) Wilcoxon test  
% h = 0, fail to reject the null hypothesis, h = 1 reject null
%N2HRC= concantenated replicates 1 and 2 hermaphrodite

%For N2 hermpahrodite replicates overlap
if exist('perATSoverlap_N2HR1', 'var') && exist('perATSoverlap_N2HR2', 'var')
    [h_N2HRATS_overlap, p_N2HRATS_overlap]= ttest2(perATSoverlap_N2HR1, perATSoverlap_N2HR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)

    N2HRC_ATSoverlap = vertcat(perATSoverlap_N2HR1, perATSoverlap_N2HR2);
end

%For N2 hermpahrodite replicates intron only
if exist('perATSIntron_N2HR1', 'var') && exist('perATSIntron_N2HR2', 'var')
    [h_N2HRATS_intron, p_N2HRATS_intron]= ttest2(perATSIntron_N2HR1, perATSIntron_N2HR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)

    N2HRC_ATSintron = vertcat(perATSIntron_N2HR1, perATSIntron_N2HR2);
end



%For dual tagged hermaphrodite replicates
if exist('perATSoverlap_DTHR1', 'var') && exist('perATSoverlap_DTHR2', 'var')
    [h_DTATS_overlap, p_DTATS_overlap]= ttest2(perATSoverlap_DTHR1, perATSoverlap_DTHR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)

    DTHRC_ATSoverlap = vertcat(perATSoverlap_DTHR1, perATSoverlap_DTHR2);
end

%For dual tagged hermpahrodite replicates intron only
if exist('perATSIntron_DTHR1', 'var') && exist('perATSIntron_DTHR2', 'var')
    [h_DTHRATS_intron, p_DTHRATS_intron]= ttest2(perATSIntron_DTHR1, perATSIntron_DTHR2);
    %[p_N2HRATSW, h_N2HRATSW] = ranksum(perATSpCR_N2HR1, perATSpCR_N2HR2)

    DTHRC_ATSintron = vertcat(perATSIntron_DTHR1, perATSIntron_DTHR2);
end




%% Comparison of N2 hermaphrodites and Dual tagged hermaphrodites

%for overlap
if exist('N2HRC_ATSoverlap', 'var') && exist('DTHRC_ATSoverlap', 'var')
    [h_N2vsDTH_overlap, p_N2vsDTH_overlap] = ttest2(N2HRC_ATSoverlap,DTHRC_ATSoverlap);
    mean_N2HC_overlap = mean(N2HRC_ATSoverlap);
    std_N2HC_overlap = std(N2HRC_ATSoverlap);
    mean_DTHC_overlap = mean(DTHRC_ATSoverlap);
    std_DTHC_overlap = std(DTHRC_ATSoverlap);
    %Graph of overlap only
    figure
    hold on
    hbar = bar(xcr, [mean_N2HC_overlap',mean_DTHC_overlap']);
    pause(0.1)
    errorbar(xcr+hbar(1).XOffset,mean_N2HC_overlap',std_N2HC_overlap'/sqrt(25), 'k.', 'linewidth', 2);
    errorbar(xcr+hbar(2).XOffset,mean_DTHC_overlap',std_DTHC_overlap'/sqrt(25), 'k.', 'linewidth', 2);
    legend('N2 overlap' ,'dual tagged mpk-1 overlap')
    hold off
end
    

%for intron only

if exist('N2HRC_ATSintron', 'var') && exist('DTHRC_ATSintron', 'var')
    [h_N2vsDTH_intron, p_N2vsDTH_intron] = ttest2(N2HRC_ATSintron,DTHRC_ATSintron);
    mean_N2HC_intron = mean(N2HRC_ATSintron);
    std_N2HC_intron = std(N2HRC_ATSintron);
    mean_DTHC_intron = mean(DTHRC_ATSintron);
    std_DTHC_intron = std(DTHRC_ATSintron);
    %Graph of intron only
    figure
    hold on
    hbar = bar(xcr, [mean_N2HC_intron',mean_DTHC_intron']);
    pause(0.1)
    errorbar(xcr+hbar(1).XOffset,mean_N2HC_intron',std_N2HC_intron'/sqrt(25), 'k.', 'linewidth', 2);
    errorbar(xcr+hbar(2).XOffset,mean_DTHC_intron',std_DTHC_intron'/sqrt(25), 'k.', 'linewidth', 2);
    legend('N2 intron only', 'dual tagged mpk-1 intron only')
    hold off
end

%dual tagged intron vs overlap
figure
hold on
hbar = bar(xcr,[mean_DTHC_overlap',mean_N2HC_intron']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,mean_DTHC_overlap',std_DTHC_overlap'/sqrt(25), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,mean_N2HC_intron',std_DTHC_intron'/sqrt(25), 'k.', 'linewidth', 2);
legend('Intron exon overlap', 'intron only')
hold off

%N2 intron vs exon
figure
hold on
hbar = bar(xcr,[mean_N2HC_overlap',mean_N2HC_intron']);
pause(0.1)
errorbar(xcr+hbar(1).XOffset,mean_N2HC_overlap',std_N2HC_overlap'/sqrt(50), 'k.', 'linewidth', 2);
errorbar(xcr+hbar(2).XOffset,mean_N2HC_intron',std_N2HC_intron'/sqrt(50), 'k.', 'linewidth', 2);
legend('Intron exon overlap', 'intron only')
hold off

    %Graph of intron/exon everything
    %has issues getting all error bars on.
    figure
    hold on
    hbar = bar(xcr, [mean_N2HC_overlap',mean_DTHC_overlap',mean_N2HC_intron',mean_DTHC_intron']);
    pause(0.1)
    errorbar(xcr+hbar(1).XOffset,mean_N2HC_overlap',std_N2HC_overlap'/sqrt(25), 'k.', 'linewidth', 2);
    pause(0.1)
    errorbar(xcr+hbar(2).XOffset,mean_DTHC_overlap',std_DTHC_overlap'/sqrt(25), 'k.', 'linewidth', 2);
    pause(0.1)
    errorbar(xcr+hbar(3).XOffset,mean_N2HC_intron',std_N2HC_intron'/sqrt(25), 'k.', 'linewidth', 2);
    pause(0.1)
    errorbar(xcr+hbar(4).XOffset,mean_DTHC_intron',std_DTHC_intron'/sqrt(25), 'k.', 'linewidth', 2);
    legend('N2 overlap', 'Dual tagged overlap', 'N2 intron only', 'Dual tagged intron only', 'location', 'northwest')
    hold off






%% Comparison between genotypes for ATS in space
% 2 tailed unpaired t-test
%put in lines to make figure

%N2 vs dual tagged hermahrodites
if exist('N2HRC_ATSpercent', 'var') && exist('DTHRC_ATSpercent', 'var') 
    [h_N2vsDTH, p_N2vsDTH]= ttest2(N2HRC_ATSpercent, DTHRC_ATSpercent);
    meanN2HRC_ATSpercent = mean(N2HRC_ATSpercent);
    meanDTHRC_ATSpercent = mean(DTHRC_ATSpercent);
    stdN2HRC_ATSpercent = std(N2HRC_ATSpercent);
    stdDTHRC_ATSpercent = std(DTHRC_ATSpercent);
    [h_N2vsDTH_ATSpercent, p_N2vsDTH_ATSpercent]=ttest2(meanN2HRC_ATSpercent,meanDTHRC_ATSpercent);
    figure
    hold on
    hbar=bar(xcr,[meanN2HRC_ATSpercent',meanDTHRC_ATSpercent']);
    pause(0.1)
    errorbar(xcr+hbar(1).XOffset,meanN2HRC_ATSpercent',stdN2HRC_ATSpercent'/sqrt(25), 'k.', 'linewidth', 2);
    errorbar(xcr+hbar(2).XOffset,meanDTHRC_ATSpercent',stdDTHRC_ATSpercent'/sqrt(25), 'k.', 'linewidth', 2);
    legend('N2 hermaprodite', 'Dual tagged mpk-1 hermaphrodite', 'location', 'northeast');
    xticks(0:10:100)
    hold off


end


%% boxplot for # mRNA of NTS+ or NTS- nuclei (modify for I/E overlap)
% 
% ntsPosNuc = cell(lic,1);
% ntsNegNuc = cell(lic,1);
% 
% for i=1:lic
%     ntsPosNuc{i} = aff{i,1}(aff{i,3}(:,8) > 0,:);
%     ntsNegNuc{i} = aff{i,1}(aff{i,3}(:,8) == 0,:);
% end
% 
% ntsPosNuc = cell2mat(ntsPosNuc);
% ntsNegNuc = cell2mat(ntsNegNuc);
% 
% ntsPosNuc(ntsPosNuc(:,1) > 20/iminfo(6),:) =[];
% ntsNegNuc(ntsNegNuc(:,1) > 20/iminfo(6),:) =[];
% 
% 
% 
% figure
% met =     11    ;
% bxplot(ntsPosNuc(:,met), ntsNegNuc(:,met));
% ylabel('# mRNA in ROI (2.5 \mu\itm\rm)', 'fontsize',15, 'interpreter', 'tex');
% set(gca,'XTickLabel',{'NTS+', 'NTS-'})
% 
% [a,b]=ttest2(ntsPosNuc(:,met), ntsNegNuc(:,met))
%% ============== # ATS per cell in space  =============================
%file name, ATS channel, y scale, iminfo
%fix axis titles

NumbATSperCell(aff, 1, 4, iminfo)

NumbATSperCell(aff_DTHR1, 1, 4, iminfo_DTHR1)

%% Intron/exon overlap



%%%--------------- categorize nuc by #NTS
anuc1 = atsNucAll(atsNucAll(:,8) == 1 & atsNucAll(:,1) < 30 / iminfo_N2HR1(6),:);
anuc2 = atsNucAll(atsNucAll(:,8) == 2 & atsNucAll(:,1) < 30 / iminfo_N2HR1(6),:);
anuc3 = atsNucAll(atsNucAll(:,8) == 3 & atsNucAll(:,1) < 30 / iminfo_N2HR1(6),:);
anuc4 = atsNucAll(atsNucAll(:,8) == 4 & atsNucAll(:,1) < 30 / iminfo_N2HR1(6),:);
anuc0 = atsNucAll(atsNucAll(:,8) == 0 & atsNucAll(:,1) < 30 / iminfo_N2HR1(6),:);
anupl = atsNucAll(atsNucAll(:,8) > 0 & atsNucAll(:,1) < 30 / iminfo_N2HR1(6),:);

%%% count # nuc with certain # of NTS.
NTSnumPerNuc = zeros(lic,6);

for i=1:lic
    NTSnumPerNuc(i,:) = [length(anuc0(anuc0(:,10) == i,10)) length(anuc1(anuc1(:,10) == i,10)) ...
        length(anuc2(anuc2(:,10) == i,10)) length(anuc3(anuc3(:,10) == i,10)) ...
        length(anuc4(anuc4(:,10) == i,10)) length(anupl(anupl(:,10) == i,10))];
end

PerNuc = NTSnumPerNuc;    % NTSnumPerNuc in percentage
totN = PerNuc(:,1) + PerNuc(:,6);
PerNuc(:,1) = PerNuc(:,1) ./ totN * 100;
PerNuc(:,2) = PerNuc(:,2) ./ totN * 100;
PerNuc(:,3) = PerNuc(:,3) ./ totN * 100;
PerNuc(:,4) = PerNuc(:,4) ./ totN * 100;
PerNuc(:,5) = PerNuc(:,5) ./ totN * 100;
PerNuc(:,6) = PerNuc(:,6) ./ totN * 100;
