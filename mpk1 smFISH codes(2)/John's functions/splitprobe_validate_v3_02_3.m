function [prna, nuc] = splitprobe_validate_v3_02_3(prna, nuc, ChOrder, iminfo, ...
    c1, c1_isIntron, c1_outcolumn, c2, c2_isIntron, c2_outcolumn)
    %% Set parameters
    distATS =       1.5            ;   % acceptible distance (um) of ATSs on different channels.
    %distATS was originally 1.5, found 18/10/03 that this is optimal
    
    %% Assert that c1 and c2 are valid channel numbers
    assert(ismember(c2, ChOrder) && ismember(c1, ChOrder), "Trying to validate invalid channels")
    
    %% set up a struct to hold channel information
    c = struct('id',{},'isIntron', {},'outcolumn',{}, 'data', {},...
        'mean_intensity', {}); 
    c(1) = setfield(c,'id',c1,'isIntron', c1_isIntron,'outcolumn',c1_outcolumn, 'prna', prna{ChOrder == c1});
    c(2) = setfield(c,'id',c2,'isIntron', c2_isIntron,'outcolumn',c2_outcolumn, prna', prna{ChOrder == c2});
    
    %% If either prna doesnt exist, return because validation is futile
    if isempty(c(1).prna) || isempty(c(2).prna)
        return
    end
    %% remove cytoplasmic mRNAs and create mean intensity values
    for n = 1:2
        if ~c(n).isIntron && ~isempty(c(n).prna)
            TF1 = (c(n).prna(:,4)==0);%TF1 = nuclear ID =0 (or cytoplasmic mRNA spot)
            c(n).prna(TF1,:) = [];
        end
        if ~isempty(c(n).prna)
            c(n).mean_intensity = mean(c(n).prna(:,6));
        else
            c(n).mean_intensity = 0;
        end
    end
    %% bulk of the analysis
    m = 0; % m is a stand-in for the opposite channel
    for n = 1:2
        if n == 1; m = 2; else; m = 1; end;
        for i = 1:length(c(n).prna(:,1)) %take the length of x coordinates
            zpln = c(n).prna(i,7);%zpln is z plane information from prna
            m_prnas = c(m).prna(c(m).prna(:,7) > zpln- 6 & c(m).prna(:,7) < zpln+ 6,:); %DetRext only for detected signals from exon channel
            %z-plane > z-plane intron-6 and z-plane of exon but < z-plane of intron
            distMat = zeros(length(m_prnas(:,1)),6);%make distMat with length of DetRext2 and 5 columns
            distMat(:,1) = m_prnas(:,1) - c(n).prna(i,1);% x coordinate intron, first column
            distMat(:,2) = m_prnas(:,2) - c(n).prna(i,2);%y coordinate intron, second column
            distMat(:,3) = m_prnas(:,6);%total signal intensity from exon
            distMat(:,6) = m_prnas(:,7) - zpln1; %distMat21(:,6) is now the z plane difference
            %                 distMat(:,3) = DetRext(:,3) - c(n).prna(i,3);
            distMat(:,4) = sqrt((distMat(:,1)*iminfo(5)).^2 + (distMat(:,2)*iminfo(6)).^2);% + (distMat(:,6)*iminfo(7)).^2); % + distMat(:,3).^2); %hypotenuse of x and y and z
            %validate 3 distances
            distMat(:,5) = m_prnas(:,8);%RNA ID from exon
            distMats = distMat(distMat(:,4) <= distATS, :);%distMats hypotenuse must be <= previously defined acceptable ATS distance
            if isempty(distMats)%if datapoint is absent from distMats, but present in intron dataset (based on prevous line)
                c(n).prna(i,c(n).outcolumn) = c(n).id*111; %no match between intron and exon
                continue
            end
            distMats = sortrows(distMats,-3);%sort least to greatest by signal intensity %JOHN: did not test

            %Test candidates selected from distMats for matching
            if distMats(1,4) <= distATS%if hypotenuse of distMats <= previously defined acceptable ATS distance
                if distMats(1,3) < mRNAIntM/2 && c(n).prna(i,6) < pRNAIntM/2 && c(n).prna(i,1) > min(c(n).prna(:,1)) + 25/iminfo(6)
                    disp(distMats(1,3))
                    %Disqualifying conditions for match:
                    %1) exon signal < mean exon signal/2
                    %2) intron signal < mean intron signal/2
                    %3) x distance intron > min x distance intron+25/iminfo(6)
                    c(n).prna(i,c(n).outcolumn) = c(n).id*111;  % mark intron only.
                    countp = countp + 1; %don't know what this is for
                else % mark ATS as matching the exon channel
                    %c(n).prna(ExATSrna_ATS(:,8) == distMats(1,5),9) = 1212;%prna{exon} of RNA ID equals RNA ID in distMats, put match in 9th column of intron
                    %%%                 table exon RNA_ID    detected exon RNA_ID

                    c(n).prna(i,c(n).outcolumn) = c(n).id*1010+c(m).id*101; 
                    c(n).prna(i,c(n).outcolumn+1) = distMats(1,5);
                    c(n).prna(i,c(n).outcolumn+2) = distMats(1,3);%put exon singnal intensity into 10th column of intron prna{ChOrder==1}
                    countm = countm + 1;%don't know what this is for
                end
            else
                c(n).prna(i,c(n).outcolumn) = -999; %if no matches with previous lines, not ATS
            end
        end
    end
    prna{ChOrder == c2} = c(2).prna;
    prna{ChOrder == c2} = c(2).prna;
    
end
