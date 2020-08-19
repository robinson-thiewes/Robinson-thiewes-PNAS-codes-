    %% cross-validate trx sites by Intron probe and by Exon probe.
    %from smFISH_v1_02_mpk1_v1 by Sarah and Changhwan, editted by John
function [ExATSrna_ATS, prna, ExATSrna] = validate(prna, ExATSrna, ChOrder, iminfo, nuc)

    % Also remove mRNA spots that are determined as ATS.
    
    % 'DexRext' 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
    %               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
    %               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
    %               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|
    
    if ismember(2, ChOrder)
        if ~isempty(ExATSrna{ChOrder == 2})%if ExATSrna Channel 2 not empty %John: channel 2 is exon
            ExATSrna_ATS= ExATSrna{ChOrder == 2}; %rename array and make matrix
            TF1 = (ExATSrna{ChOrder == 2}(:,4)==0); %TF1 = nuclear ID =0 (or cytoplasmic mRNA spot)
            ExATSrna_ATS(TF1,:) = []; %remove cytoplasmic detected RNAs
            mRNAIntM = mean(ExATSrna_ATS(:,6)); % mean intensity of all exon detected ATS spots.
        else
            mRNAIntM = 0;
        end
    else
        mRNAIntM = 0;
    end
    
    if ismember(1, ChOrder)
        if ~isempty(prna{ChOrder == 1})
            %SARAH- if prna has values
            % associated with a nucleus that we need to worry about? SARAH-
            % already checked for nuclear association
            pRNAIntM = mean(prna{ChOrder == 1}(:,6)); %mean intensity of intron spots
        else
            pRNAIntM = 0;
        end
    else
        pRNAIntM = 0;
    end
    
    if ismember(2, ChOrder) && ismember(1, ChOrder)
        prna{ChOrder == 1}(:,8:11) = 0;
        ExATSrna_ATS(:,10:12) = 0;
        distATS =       0.5            ;   % acceptible distance (um) of ATSs on different channels.
        
        countp = 0;
        countm = 0;
        
        % remove detected ATS spots from the list of mRNAs
        if ~isempty(prna{ChOrder == 1})%if prna{ChOrder == 1} is not empty
            for i = 1:length(prna{ChOrder == 1}(:,1)) %take the length of x coordinates
                zpln1 = prna{ChOrder == 1}(i,7);%zpln is z plane information from prna
                DetRext2 = ExATSrna_ATS(ExATSrna_ATS(:,7) > zpln1- 6 & ExATSrna_ATS(:,7) < zpln1+ 6,:); %DetRext only for detected signals from exon channel
                %z-plane > z-plane intron-6 and z-plane of exon but < z-plane of intron
                distMat12 = zeros(length(DetRext2(:,1)),5);%make distMat12 with length of DetRext2 and 5 columns
                distMat12(:,1) = DetRext2(:,1) - prna{ChOrder == 1}(i,1);% x coordinate intron, first column
                distMat12(:,2) = DetRext2(:,2) - prna{ChOrder == 1}(i,2);%y coordinate intron, second column
                distMat12(:,3) = DetRext2(:,6);%total signal intensity from exon
                %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder == 1}(i,3);
                distMat12(:,4) = sqrt(distMat12(:,1).^2 + distMat12(:,2).^2); % + distMat(:,3).^2); %hypotenuse of x and y
                distMat12(:,5) = DetRext2(:,8);%RNA ID from exon
                distMats12 = distMat12(distMat12(:,4) <= distATS/iminfo(6), :);%distMats12 hypotenuse must be <= previously defined acceptable ATS distance
                if isempty(distMats12)%if datapoint is absent from distMats12, but present in intron dataset (based on prevous line)
                    prna{ChOrder == 1}(i,8) = 111; %no match between intron and exon
                    continue
                end
                distMats12 = sortrows(distMats12,-3);%sort least to greatest by signal intensity %JOHN: did not test
                
                %Test candidates selected from distMats for matching
                if distMats12(1,4) <= distATS/iminfo(6)%if hypotenuse of distMats12 <= previously defined acceptable ATS distance
                    if distMats12(1,3) < mRNAIntM/2 && prna{ChOrder == 1}(i,6) < pRNAIntM/2 && prna{ChOrder == 1}(i,1) > min(prna{ChOrder == 1}(:,1)) + 25/iminfo(6)
                        %Disqualifying conditions for match:
                        %1) exon signal < mean exon signal/2
                        %2) intron signal < mean intron signal/2
                        %3) x distance intron > min x distance intron+25/iminfo(6)
                        prna{ChOrder == 1}(i,8) = 111;  % mark intron only.
                        countp = countp + 1; %don't know what this is for
                    else % mark ATS as matching the exon channel
                        %prna{ChOrder == 1}(ExATSrna_ATS(:,8) == distMats12(1,5),9) = 1212;%prna{exon} of RNA ID equals RNA ID in distMats12, put match in 9th column of intron
                        %%%                 table exon RNA_ID    detected exon RNA_ID
                        prna{ChOrder == 1}(i,9) = 1212;
                        prna{ChOrder == 1}(i,10) = distMats12(1,3);%put exon singnal intensity into 10th column of intron prna{ChOrder==1}
                        countm = countm + 1;%don't know what this is for
                    end
                else
                    prna{ChOrder == 1}(i,11) = 999; %if no matches with previous lines, not ATS
                end
            end
        end
        
        %exon vs intron comparison
        if ~isempty(ExATSrna_ATS) %for exon channel ATS
            for i = 1:length(ExATSrna_ATS(:,1)) %i=length of x-coordinates
                zpln2 = ExATSrna_ATS(i,7); %zpln2 is prna(exonATS) z coordinate in plane #
                DetRext1 = prna{ChOrder == 1}(prna{ChOrder == 1}(:,7) > zpln2- 6 & prna{ChOrder == 1}(:,7) < zpln2+ 6,:);
                % DetRext1 = prna(1st) z coordinate > 2nd half z-coordinate - 6 logical AND 1st z-coordinate < 2nd half z-coordinate + 6 for all columns
                distMat21 = zeros(length(DetRext1(:,1)),5); % set up distMat21 array to be zeros for the length of DetRext1(all, x-coordinate), 5 columns
                distMat21(:,1) = DetRext1(:,1) - ExATSrna_ATS(i,1); %distMat21(exon vs intron)= intron x-coordinate - prna(exon) x-coordinate
                distMat21(:,2) = DetRext1(:,2) - ExATSrna_ATS(i,2); %distMat21(exon vs intron)= intron y-coordinate - prna(exon) y-coordinate
                distMat21(:,3) = DetRext1(:,6); %distMat21(exon vs intron)= intron signal intensity
                %                 distMat(:,3) = DetRext(:,3) - prna{ChOrder ==1}(i,3); This is for z-coordinate, not used
                distMat21(:,4) = sqrt(distMat21(:,1).^2 + distMat21(:,2).^2); % + distMat21(:,3).^2); distMat21 (exon vs intron) = square root of x^2 and y^2 coordinate distMat21
                distMat21(:,5) = DetRext1(:,4); % distMat21(exon vs intron) = intron nuclear #
                distMats21 = distMat21(distMat21(:,4) <= distATS/iminfo(6), :); %distMat21(exon vs intron) = distMat21(exon vs intron) hypotenuse <= distanceATS/iminfo6 (um measuremet)
                
                if isempty(distMats21) %if distMats21(exon vs intron) empty
                    ExATSrna_ATS(i,10) = 222; %exon probe and intron probe ATS do not match
                    continue
                end
                
                distMats21 = sortrows(distMats21,-3); %sort distMat21(exon vs exon) based on signal intensity, least to greatest
                
                %Test candidates selected from distMats for matching
                if distMats21(1,4) <= distATS/iminfo(6)%if hypotenuse of distMats21 <= previously defined acceptable ATS distance
                    if distMats21(1,3) < pRNAIntM/2 && ExATSrna_ATS(i,6) < mRNAIntM/2 && ExATSrna_ATS(i,1) > min(ExATSrna_ATS(:,1)) + 25/iminfo(6)
                        %Disqualifying conditions for match:
                        %1) intron signal < mean intron signal/2
                        %2) exon signal < mean exon signal/2
                        %3) x distance exon > min x distance exon+25/iminfo(6)
                        ExATSrna_ATS(i,10) = 222;  % mark exon spot only.
                        countp = countp + 1; %don't know what this is for
                    else % mark ATS in intron channel
                        %ExATSrna_ATS(prna{ChOrder == 1}(:,4) == distMats21(1,5),9) = 2121;%prna{intron} nuclear ID equals nuclear ID in distMats12, put match in 13th column of exon
                        ExATSrna_ATS(i,9) = 2121;
                        ExATSrna_ATS(i,11) = distMats21(1,3);%put matching intron intensity into exon prna{ChOrder==2}
                        countm = countm + 1;%don't know what this is for
                    end
                else
                    ExATSrna_ATS(i,12) = 999; %if no matches with previous lines, not ATS
                end
            end
        end
    end
    % Comment out if intron/exon ratio needs to be analyzed
    %             % remove intron ATS that are not seen at exon channel
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
    
    
    %uncomment out if intron channel
    %%%%%     Record ATS intensity    %%%%%
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
end