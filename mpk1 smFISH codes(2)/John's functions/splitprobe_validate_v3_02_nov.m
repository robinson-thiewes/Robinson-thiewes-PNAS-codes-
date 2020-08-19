%% cross-validate trx sites by Intron probe and by Exon probe.
function [prna, nuc] = splitprobe_validate_v3_02(prna, nuc, ChOrder, iminfo)


% Also remove mRNA spots that are determined as ATS.

% 'DexRext' 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
%               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
%               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
%               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|

%for exons take out because variables not needed
%         if ismember(2, ChOrder)
%             if ~isempty(prna{ChOrder == 2})
%                 mRNAIntM = mean(prna{ChOrder == 2}(:,6)); % mean intensity of all mRNA spots.
%             else
%                 mRNAIntM = 0;
%             end
%         else
%             mRNAIntM = 0;
%         end
%
%         %For intron first half
%         if ismember(1, ChOrder)
%             if ~isempty(prna{ChOrder == 1})
%                 pRNAIntM1 = mean(prna{ChOrder == 1}(:,6));
%             else
%                 pRNAIntM1 = 0;
%             end
%         else
%             pRNAIntM1 = 0;
%         end
%
%         %For intron second half
%         if ismember(3, ChOrder)
%             if ~isempty(prna{ChOrder == 3})
%                 pRNAIntM2 = mean(prna{ChOrder == 3}(:,6));
%             else
%                 pRNAIntM2 = 0;
%             end
%         else
%             pRNAIntM2 = 0;
%         end

%cross validate for all three channels

if ismember(2, ChOrder) && ismember(1, ChOrder) && ismember(3, ChOrder)
    prna{ChOrder == 1}(:,8:13) = 0; %1st half
    prna{ChOrder == 2}(:,13:18) = 0;%exon
    prna{ChOrder == 3}(:,8:13) = 0;%2nd half
    distATS =       0.5            ;   % acceptible distance (um) of ATSs on different channels
    %             countp = 0; %comment out because not using part of code that removes ATS not seen in exon channel
    %             countm = 0;
    
end


% remove detected ATS spots from the list of mRNAs
% for 1st half vs exon

if ~isempty(prna{ChOrder == 1}) %if prna(intron 1st half) is NOT empty
    %This loop compares 1st and exon channels
    for i = 1:length(prna{ChOrder == 1}(:,1)) %i=length of x-coordinates
        zpln1 = prna{ChOrder == 1}(i,7); %zpln1 is prna(1st half) z coordinate in plane #
        DetRext2 = prna{ChOrder == 2}(prna{ChOrder == 2}(:,7) > zpln1- 6 & prna{ChOrder == 2}(:,7) < zpln1+ 6,:);
        % DetRext2 = prna(exon) z coordinate > 1st half z-coordinate - 6 logical AND exon z-coordinate < 1st half z-coordinate + 6 for all columns
        distMat12 = zeros(length(DetRext2(:,1)),5); % set up distMat12 array to be zeros for the length of DetRext2(all, x-coordinate), 5 columns
        distMat12(:,1) = DetRext2(:,1) - prna{ChOrder == 1}(i,1); %distMat12 (1st vs exon)= 1st half x-coordinate - prna(2nd) x-coordinate
        distMat12(:,2) = DetRext2(:,2) - prna{ChOrder == 1}(i,2); %distMat12(1st vs exon)= 1st half y-coordinate - prna(2nd) y-coordinate
        distMat12(:,3) = DetRext2(:,6); %distMat12(1st vs exon)= total signal intensity for 1st half
        %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder ==1}(i,3); This is for z-coordinate, not used
        distMat12(:,4) = sqrt(distMat12(:,1).^2 + distMat12(:,2).^2); % + distMat12(:,3).^2); distMat12 (1st vs exon) = square root of x^2 and y^2 coordinate distMat1
        distMat12(:,5) = DetRext2(:,8); % distMat2(1st vs exon) = RNA ID
        distMats12 = distMat12(distMat12(:,4) <= distATS/iminfo(6), :); %distMat12(1st vs exon) = distMat12(1st vs exon) hypotenuse <= distanceATS/iminfo6 (um measuremet)
        
        if isempty(distMats12) %if distMats12(1st vs exon) empty
            prna{ChOrder == 1}(i,9) = 111; % 1st vs exon ATS do not match
            continue
        end
        
        if ~isempty(ExATSrna{ChOrder == 2})%if ExATSrna Channel 2 not empty %John: channel 2 is exon
            distMats12 = sortrows(distMats12,-3); %sort distMat1(1st vs exon) based on signal intensity, least to greatest
        end
        if distMats12(1,4) <= distATS/iminfo(6)%if hypotenuse is <= distance ATS/iminfo (um)
            prna{ChOrder == 2}(prna{ChOrder == 2}(:,4) == distMats12(1,5),9) = 1212;%If exon (all exon spots, RNA ID) equals (==) distMats12(1st vs exon)RNA ID,
            % mark ATS in exon with 1212 matches 1st half
            prna{ChOrder == 1}(i,8) = distMats12(1,3); %prna(1st half) total signal intensity = distMats(1st vs exon) total signal intensity
            prna{ChOrder == 1}(i,10) = distMats12(1,5); %put RNA ID into column 10
            %                             countm = countm + 1; %see above
            
        else
            prna{ChOrder == 1}(i,9) = 1212; %put exon matching into prna(1st) into column 9
        end
    end
    %This loop compares 1st vs 2nd half channels
    for i = 1:length(prna{ChOrder == 1}(:,1)) %i=length of x-coordinates
        zpln1 = prna{ChOrder == 1}(i,7); %zpln1 is prna(1st half) z coordinate in plane #
        DetRext3 = prna{ChOrder == 3}(prna{ChOrder == 3}(:,7) > zpln1- 6 & prna{ChOrder == 3}(:,7) < zpln1+ 6,:);
        % DetRext3 = prna(2nd) z coordinate > 1st half z-coordinate - 6 logical AND 2nd z-coordinate < 1st half z-coordinate + 6 for all columns
        distMat13 = zeros(length(DetRext3(:,1)),5); % set up distMat13 array to be zeros for the length of DetRext3(all, x-coordinate), 5 columns
        distMat13(:,1) = DetRext3(:,1) - prna{ChOrder == 1}(i,1); %distMat13(1st vs 2nd)= 2nd x-coordinate - prna(2nd half) x-coordinate
        distMat13(:,2) = DetRext3(:,2) - prna{ChOrder == 1}(i,2); %distMat13(1st vs 2nd)= 2nd y-coordinate - prna(2nd half) y-coordinate
        distMat13(:,3) = DetRext3(:,6); %distMat13(1st vs 2nd)= 2nd signal intensity
        %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder ==1}(i,3); This is for z-coordinate, not used
        distMat13(:,4) = sqrt(distMat13(:,1).^2 + distMat13(:,2).^2); % + distMat13(:,3).^2); distMat13 (1st vs 2nd) = square root of x^2 and y^2 coordinate distMat13
        distMat13(:,5) = DetRext3(:,4); % distMat13(1st vs 2nd) = nuclear number
        distMats13 = distMat13(distMat13(:,4) <= distATS/iminfo(6), :); %distMat13(1st vs 2nd) = distMat13(1st vs 2nd) hypotenuse <= distanceATS/iminfo6 (um measuremet)
        
        if isempty(distMats13) %if distMats13(1st vs 2nd) empty
            prna{ChOrder == 1}(i,12) = 111; % 2nd probe and 1st probe ATS do not match
            continue
        end
        
        distMats13 = sortrows(distMats13,-3); %sort distMat13(1st vs 2nd) based on signal intensity, least to greatest
        
        if distMats13(1,4) <= distATS/iminfo(6)
            prna{ChOrder == 3}(prna{ChOrder == 3}(:,8) == distMats13(1,5),19) = 1313;%If 2nd (all 2nd spots, nuc #) equals (==) distMats13(1st vs 2nd) nuclear number,
            % mark ATS in 2nd channel with 1313 matches 2nd half
            prna{ChOrder == 1}(i,11) = distMats13(1,3); %prna(1st half) total signal intensity = distMats(1st vs 2nd) total signal intensity
            prna{ChOrder == 1}(i,13) = distMats13(1,5);%Put nuclear # into 1st half column 13
            %countm = countm + 1; %see above
            
        else
            prna{ChOrder == 1}(i,12) = 1313; %put 2nd matching into prna(2nd) into column 12
        end
        
    end
end


%for 2nd half probe
if ~isempty(prna{ChOrder == 3}) %if 2nd half not empty
    %This loop compares 2nd vs exon
    for i = 1:length(prna{ChOrder == 3}(:,1)) %i=length of x-coordinates
        zpln3 = prna{ChOrder == 3}(i,7); %zpln3 is prna(2nd half) z coordinate in plane #
        DetRext2 = prna{ChOrder == 2}(prna{ChOrder == 2}(:,7) > zpln3- 6 & prna{ChOrder == 2}(:,7) < zpln3+ 6,:);
        % DetRext2 = prna(exon) z coordinate > 2nd half z-coordinate - 6 logical AND exon z-coordinate < 2nd half z-coordinate + 6 for all columns
        distMat32 = zeros(length(DetRext2(:,1)),5); % set up distMat32 array to be zeros for the length of DetRext2(all, x-coordinate), 5 columns
        distMat32(:,1) = DetRext2(:,1) - prna{ChOrder == 3}(i,1); %distMat32 (2nd vs exon)= exon x-coordinate - prna(2nd) x-coordinate
        distMat32(:,2) = DetRext2(:,2) - prna{ChOrder == 3}(i,2); %distMat32(2nd vs exon)= exon y-coordinate - prna(2nd) y-coordinate
        distMat32(:,3) = DetRext2(:,6); %distMat32(2nd vs exon)= exon signal intensity
        %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder ==1}(i,3); This is for z-coordinate, not used
        distMat32(:,4) = sqrt(distMat32(:,1).^2 + distMat32(:,2).^2); % + distMat32(:,3).^2); distMat32 (2nd vs exon) = square root of x^2 and y^2 coordinate distMat23
        distMat32(:,5) = DetRext2(:,8); % distMat3(2nd vs exon) = exon RNA ID
        distMats32 = distMat32(distMat32(:,4) <= distATS/iminfo(6), :); %distMat32(2nd vs exon) = distMat32(2nd vs exon) hypotenuse <= distanceATS/iminfo6 (um measuremet)
        
        if isempty(distMats32) %if distMats23(2nd vs exon) empty
            prna{ChOrder == 3}(i,9) = 222; % 2nd half probe and exon probe ATS do not match
            continue
        end
        
        distMats32 = sortrows(distMats32,-3); %sort distMat23(2nd vs exon) based on signal intensity, least to greatest
        
        if distMats32(1,4) <= distATS/iminfo(6)
            prna{ChOrder == 2}(prna{ChOrder == 3}(:,8) == distMats32(1,5),9) = 2323;%If exon (all exon spots, RNA ID) equals (==) distMats32(2nd vs exon) RNA ID,
            % mark ATS in exon channel with 2323 matches 2nd half
            prna{ChOrder == 3}(i,8) = distMats32(1,3); %prna(exon) total signal intensity = distMats32(2nd vs exon) total signal intensity
            %                             countm = countm + 1; %see above
            prna{ChOrder == 3}(i,10) = distMats32(1,5);%Put exon RNA ID into 2nd half column 10
            
        else
            prna{ChOrder == 3}(i,9) = 2323; %put exon matching into prna(2nd) into column 9
        end
    end
    %This loop compares 2nd vs 1st intron channels
    for i = 1:length(prna{ChOrder == 3}(:,1)) %i=length of x-coordinates
        zpln3 = prna{ChOrder == 3}(i,7); %zpln3 is prna(2nd half) z coordinate in plane #
        DetRext1 = prna{ChOrder == 1}(prna{ChOrder == 1}(:,7) > zpln3- 6 & prna{ChOrder == 1}(:,7) < zpln3+ 6,:);
        % DetRext1 = prna(1st) z coordinate > 2nd half z-coordinate - 6 logical AND 1st z-coordinate < 2nd half z-coordinate + 6 for all columns
        distMat31 = zeros(length(DetRext1(:,1)),5); % set up distMat21 array to be zeros for the length of DetRext1(all, x-coordinate), 5 columns
        distMat31(:,1) = DetRext1(:,1) - prna{ChOrder == 3}(i,1); %distMat31(2nd vs 1st)= 1st x-coordinate - prna(2nd half) x-coordinate
        distMat31(:,2) = DetRext1(:,2) - prna{ChOrder == 3}(i,2); %distMat31(2nd vs 1st)= 1st y-coordinate - prna(2nd half) y-coordinate
        distMat31(:,3) = DetRext1(:,6); %distMat31(1st vs 2nd)= 2nd signal intensity
        %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder ==1}(i,3); This is for z-coordinate, not used
        distMat31(:,4) = sqrt(distMat31(:,1).^2 + distMat31(:,2).^2); % + distMat31(:,3).^2); distMat21 (2nd vs 1st) = square root of x^2 and y^2 coordinate distMat21
        distMat31(:,5) = DetRext1(:,4); % distMat31(2nd vs 1st) = 2nd nuclear #
        distMats31 = distMat31(distMat31(:,4) <= distATS/iminfo(6), :); %distMat31(2nd vs 1st) = distMat21(2nd vs 1st) hypotenuse <= distanceATS/iminfo6 (um measuremet)
        
        if isempty(distMats31) %if distMats13(2nd vs 1st) empty
            prna{ChOrder == 3}(i,12) = 222; %2nd probe and 1st probe ATS do not match
            continue
        end
        
        distMats31 = sortrows(distMats31,-3); %sort distMat13(1st vs 2nd) based on signal intensity, least to greatest
        
        if distMats31(1,4) <= distATS/iminfo(6)
            prna{ChOrder == 1}(prna{ChOrder == 1}(:,4) == distMats31(1,5),12) = 2121;%If 1st (all 1st spots, nuc #) equals (==) distMats21(2nd vs 1st) nuc #,
            % mark ATS in 2nd channel with 2121 matches 1st half in column 14
            prna{ChOrder == 3}(i,11) = distMats31(1,3); %prna(1st half) total signal intensity = distMats(1st vs exon) total signal intensity
            prna{ChOrder == 3}(i,13) = distMats31(1,5); %put nuclear ID into column 13
            %                             countm = countm + 1; %see above
            
        else
            prna{ChOrder == 3}(i,12) = 1212; %put 1st matching into prna(2nd) into column 12
        end
        
    end
end
% exon comparison
%exon vs 1st half
if ~isempty(prna{ChOrder == 2}) %for exon channel ATS
    for i = 1:length(prna{ChOrder == 2}(:,1)) %i=length of x-coordinates
        zpln2 = prna{ChOrder == 2}(i,7); %zpln2 is prna(exonATS) z coordinate in plane #
        DetRext1 = prna{ChOrder == 1}(prna{ChOrder == 1}(:,7) > zpln2- 6 & prna{ChOrder == 1}(:,7) < zpln2+ 6,:);
        % DetRext1 = prna(1st) z coordinate > 2nd half z-coordinate - 6 logical AND 1st z-coordinate < 2nd half z-coordinate + 6 for all columns
        distMat21 = zeros(length(DetRext1(:,1)),5); % set up distMat21 array to be zeros for the length of DetRext1(all, x-coordinate), 5 columns
        distMat21(:,1) = DetRext1(:,1) - prna{ChOrder == 2}(i,1); %distMat21(exon vs 1st)= 1st x-coordinate - prna(exon) x-coordinate
        distMat21(:,2) = DetRext1(:,2) - prna{ChOrder == 2}(i,2); %distMat21(exon vs 1st)= 1st y-coordinate - prna(exon) y-coordinate
        distMat21(:,3) = DetRext1(:,6); %distMat21(exon vs 1st)= 1st signal intensity
        %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder ==1}(i,3); This is for z-coordinate, not used
        distMat21(:,4) = sqrt(distMat21(:,1).^2 + distMat21(:,2).^2); % + distMat21(:,3).^2); distMat21 (exon vs 1st) = square root of x^2 and y^2 coordinate distMat21
        distMat21(:,5) = DetRext1(:,4); % distMat21(exon vs 1st) = 1st nuclear #
        distMats21 = distMat21(distMat21(:,4) <= distATS/iminfo(6), :); %distMat21(exon vs 1st) = distMat21(exon vs 1st) hypotenuse <= distanceATS/iminfo6 (um measuremet)
        
        if isempty(distMats21) %if distMats21(exon vs 1st) empty
            prna{ChOrder == 2}(i,13) = 333; %exon probe and 1st probe ATS do not match
            continue
        end
        
        distMats21 = sortrows(distMats21,-3); %sort distMat21(exon vs 1st) based on signal intensity, least to greatest
        
        if distMats21(1,4) <= distATS/iminfo(6)
            prna{ChOrder == 1}(prna{ChOrder == 1}(:,4) == distMats21(1,5),15) = 3131;%If 1st (all 1st spots, nuc #) equals (==) distMats21(exon vs 1st) nuc #,
            % mark ATS in 1st channel with 3131 matches exon in column 17
            prna{ChOrder == 2}(i,14) = distMats21(1,3); %prna(exon) total signal intensity = distMats(exon vs 1st) total signal intensity
            prna{ChOrder == 2}(i,15) = distMats21(1,5); %put nuclear ID into column 15
            %                             countm = countm + 1; %see above
            
        else
            prna{ChOrder == 2}(i,13) = 3131; %put 1st matching into prna(exon) into column 13
        end
    end
    
    %for exon vs 2nd
    for i = 1:length(prna{ChOrder == 2}(:,1)) %i=length of x-coordinates
        zpln2 = prna{ChOrder == 2}(i,7); %zpln2 is prna(exonATS) z coordinate in plane #
        DetRext3 = prna{ChOrder == 3}(prna{ChOrder == 3}(:,7) > zpln2- 6 & prna{ChOrder == 3}(:,7) < zpln2+ 6,:);
        % DetRext1 = prna(2nd) z coordinate > 2nd half z-coordinate - 6 logical AND exon z-coordinate < 2nd half z-coordinate + 6 for all columns
        distMat23 = zeros(length(DetRext3(:,1)),5); % set up distMat23 array to be zeros for the length of DetRext3(all, x-coordinate), 5 columns
        distMat23(:,1) = DetRext3(:,1) - prna{ChOrder == 2}(i,1); %distMat23(exon vs 2nd)= 2nd x-coordinate - prna(exon) x-coordinate
        distMat23(:,2) = DetRext3(:,2) - prna{ChOrder == 2}(i,2); %distMat23(exon vs 2nd)= 2nd y-coordinate - prna(exon) y-coordinate
        distMat23(:,3) = DetRext3(:,6); %distMat23(exon vs 2nd)= 2nd signal intensity
        %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder ==1}(i,3); This is for z-coordinate, not used
        distMat23(:,4) = sqrt(distMat23(:,1).^2 + distMat23(:,2).^2); % + distMat23(:,3).^2); distMat23 (exon vs 2nd) = square root of x^2 and y^2 coordinate distMat21
        distMat23(:,5) = DetRext3(:,4); % distMat23(exon vs 2nd) = 2nd nuclear #
        distMats23 = distMat23(distMat23(:,4) <= distATS/iminfo(6), :); %distMat23(exon vs 2nd) = distMat23(exon vs 2nd) hypotenuse <= distanceATS/iminfo6 (um measuremet)
        
        if isempty(distMats23) %if distMats23(exon vs 2nd) empty
            prna{ChOrder == 2}(i,16) = 333; %exon probe and 2nd probe ATS do not match
            continue
        end
        
        distMats23 = sortrows(distMats23,-3); %sort distMat23(exon vs 1st) based on signal intensity, least to greatest
        
        if distMats23(1,4) <= distATS/iminfo(6)
            prna{ChOrder == 3}(prna{ChOrder == 3}(:,4) == distMats23(1,5),16) = 3232;%If 1st (all 2nd spots, nuc #) equals (==) distMats23(exon vs 2nd) nuc #,
            % mark ATS in 2nd channel with 3131 matches exon in column 18
            prna{ChOrder == 2}(i,17) = distMats23(1,3); %prna(exon) total signal intensity = distMats(exon vs 2nd) total signal intensity
            prna{ChOrder == 2}(i,18) = distMats23(1,5); %put nuclear ID into column 15
            %                             countm = countm + 1; %see above
            
        else
            prna{ChOrder == 2}(i,16) = 3232; %put 2nd matching into prna(exon) into column 16
        end
        
    end
end









%%             % remove intron ATS that are not seen at exon channel
%commenting out because I do not want to remove ATS not
%seen in exon channel
%     %         NoOvlapATS = prna{ChOrder == 1}(prna{ChOrder == 1}(:,8) == 999,:);
%             if countp > 0
% %                 prna{ChOrder == 1}(prna{ChOrder == 1}(:,9) == 999,:) = [];
%             end
%             prna{ChOrder == 1}(:,9) = [];
%
%             if countm > 0
% %                 prna{ChOrder == 2}(prna{ChOrder == 2}(:,13) == 999,:) = [];
%             end
%             prna{ChOrder == 2}(:,13) = [];
%         end



%%%%%     Record ATS intensity first half    %%%%%
%%% 'nuc' for ATS
%       |1-3: XYZ-coordinates |4: radius (pixels)| 5: ave. circularity| 6: std circularity |
%       | 7: DAPI intensity |8: # ATS per nuc | 9: summed ATS intensity from INTRON channel |
%       | 10: summed ATS int. from EXON channel
if ~isempty(nuc{ChOrder == 1})
    nNcount = find( nuc{ChOrder == 1}(:,8) > 0 );
    for i = 1:length(nNcount)
        ints = prna{ChOrder == 1}(prna{ChOrder == 1}(:,4) == nNcount(i),6);
        nuc{ChOrder == 1}(nNcount(i),9) = sum(ints);
        ints2 = prna{ChOrder == 1}(prna{ChOrder == 1}(:,4) == nNcount(i),8);
        nuc{ChOrder == 1}(nNcount(i),10) = sum(ints2);
    end
end

%%%%%     Record ATS intensity second half    %%%%%
%%% 'nuc' for ATS
%       |1-3: XYZ-coordinates |4: radius (pixels)| 5: ave. circularity| 6: std circularity |
%       | 7: DAPI intensity |8: # ATS per nuc | 9: summed ATS intensity from INTRON channel |
%       | 10: summed ATS int. from EXON channel
if ~isempty(nuc{ChOrder == 3})
    nNcount = find( nuc{ChOrder == 3}(:,8) > 0 );
    for i = 1:length(nNcount)
        ints3 = prna{ChOrder == 3}(prna{ChOrder == 3}(:,4) == nNcount(i),6);
        nuc{ChOrder == 3}(nNcount(i),9) = sum(ints3);
        ints4 = prna{ChOrder == 3}(prna{ChOrder == 3}(:,4) == nNcount(i),8);
        nuc{ChOrder == 3}(nNcount(i),10) = sum(ints4);
    end
end
fprintf("Finished\n");
end