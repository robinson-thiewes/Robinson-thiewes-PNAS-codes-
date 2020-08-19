function [aff, mrnaNucAll, mrnaAll, atsNucAll, atsAll] = normaf(af, ChOrder, skipped) %, intron_values_all, mean_intron_values)
%%% This function aligns coordinates of all detected spots to the tip of
%%% the distal gonad. It also normalizes the RNA signal intensities based
%%% on the mean intensity of detected mRNA spots. Mean single mRNA intensity is set to 10.
%%% For input variable 'af', use 'af' from smFISH analysis.
%%% For input variable ChOrder, put   '1' for intron,     '2' for exon,   or array for multi-channel images.
%%% (e.g. [ 2    1 ] for 1st channel for Exon and 2nd channel, Intron.
%%% For input variable skipped, put number or array containing images to skip.
%%% e.g. skipped = [ 1 3], meaning 1st & 3rd images will be skipped from
%%% normalization and excluded in 'aff' output variable.
%%% Usage: [output variable] = normaf(af, [ 2  1 ], [])
%%% It outputs five variables:
%%% 1st output variable: normalized af, so called 'aff'
%%% 2nd output variable: 1st column of af, altogether in a matrix. and so
%%% forth for 3-5th variables (2nd - 4th column of af).

%%% ============= Align coordinates in 'aff' =======================

%----- METHOD 1: align by left end of the first nucleus (0 on x-axis) -----
ic = find(ChOrder == 1);
xc = find(ChOrder == 2);

aff = af;
aff(skipped,:) = [];

%x y z normalizaition, add in extra columns
lic = size(aff,1);
for i=1:lic
    if isempty(aff{i,1})
            continue
    end
    sortNuc = sort(aff{i,1}(:,1));
    minNuc = sortNuc(2);
    DTCloc = minNuc - aff{i,1}(aff{i,1}(:,1) == minNuc,4);
    if ~isempty(aff{i,1})
        aff{i,1}(:,1) = aff{i,1}(:,1) - DTCloc;
    end
    if ~isempty(aff{i,2})
        aff{i,2}(:,1) = aff{i,2}(:,1) - DTCloc;
    end
    if ~isempty(aff{i,3})
        aff{i,3}(:,1) = aff{i,3}(:,1) - DTCloc;
    end
    if ~isempty(aff{i,4})
        aff{i,4}(:,1) = aff{i,4}(:,1) - DTCloc;
    end
    if ~isempty(aff{i,5})
        aff{i,5}(:,1) = aff{i,5}(:,1) - DTCloc;
    end
    if ~isempty(aff{i,12}) %normalize for ExATSrna_ATS
        aff{i,12}(:,1) = aff{i,12}(:,1) - DTCloc;
    end
    if ~isempty(aff{i,14})%normalize ExATSnuc
        aff{i,14}(:,1) = aff{i,14}(:,1) - DTCloc;
    end
    if ~isempty(aff{i,15})%normalize ExATSrna
        aff{i,15}(:,1) = aff{i,15}(:,1) - DTCloc;
    end
end


%%%----------------- Align cyl_axis 0-end ---------------------------
for i = 1:lic
    if isempty(aff{i,5})
            continue
    end
    if aff{i,5}(1,1) ~= 0
        
        aff{i,5} = flipud(aff{i,5});
    end
end

%%%--- Normalize intensity of single spot in individual images for summary---
for i = 1:lic
    if isempty(aff{i,xc*2})
        continue
    end
    
    cytoExon_spots = aff{i,xc*2}(:,12)==0; %makes list based on distance to closest nucleus
    mRNA_spots = aff{i,xc*2}(:,:);
    mRNA_spots(cytoExon_spots,:) = []; %removes exon detection that are within nucleus
    
    foldnum = aff{i,xc*2}(:,6) / mean(mRNA_spots(:,6));
%     intron_foldnum = intron_values_all(:,6) / mean_intron_values;
    if aff{i,xc*2}(:,6) > mean(mRNA_spots(:,6)) * 4
        aff{i,xc*2}(:,6) = aff{i,xc*2}(:,6) / foldnum; 
    end
    
%     if intron_values_all(:,6) > mean_intron_values *4
%         intron_values_all(:,6) = intron_values_all(:,6)%/ intron_foldnum;
%     end
    
    % normalize NTS intensity using mean of single mRNA intensity
    if isempty(aff{i,ic*2})
        continue
    end
    %aff{i,ic*2-1}(:,10) = aff{i,ic*2-1}(:,10) / mean(aff{i,xc*2}(:,6)) *10; % from line of code that was removed for intron channel processing

%     aff{i,ic*2}(:,6) = aff{i,ic*2}(:,6) / mean_intron_values * 10;%normalize intron detections
%     aff{i,ic*2}(:,10)= aff{i,ic*2}(:,10)/mean(aff{i,xc*2}(:,6)) * 10; %normalize exon matches to single mRNA
    aff{i,xc*2}(:,6) = aff{i,xc*2}(:,6) / mean(aff{i,xc*2}(:,6)) * 10;%normalize cytoplasmic exon detections
    aff{i,12}(:,6)= aff{i,12}(:,6)/mean(aff{i,xc*2}(:,6)) *10; %normalize exon detected ATS
end


%%% label germline ID at af{i,3}
lic = size(aff,1);
for i=1:lic
    aff{i,3}(:,11) = i;
end


mrnaNucAll = cell2mat(aff(:,xc*2-1));
mrnaAll = cell2mat(aff(:,xc*2));
atsNucAll = cell2mat(aff(:,ic*2-1));
atsAll = cell2mat(aff(:,ic*2));

