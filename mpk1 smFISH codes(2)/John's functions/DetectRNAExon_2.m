function [DETrna, derna, mom] = DetectRNAExon_2(rPla, pixr, iminfo, f, lic, g, thresForExon)
% This function detects intronic RNA spots (transcription sites). The
% output is RNA coordinates, counts and intensity

% Editted from DetectRNAExon by John to include a line:
% DETrna(:,8) = DETrnaP(:,1);
% that will add RNA IDs to DETrna


sizel = iminfo(4);
fint = cell(sizel,1);  
t1 = cell(sizel,1); % cropped RNA images only in the nuclei.
mofm = zeros(sizel,1);
Pmofm = zeros(sizel,1);
pMask = cell(sizel,1);
se = strel('disk', pixr+1);
mupr = thresForExon+0.5; % multiplier to define the background level (becomes higher when bg is high)     


% peakDots = cell( sizel, 1);

% ================ %%% sharing for RNA detection ========================= 
% make general threshold
parfor i = 1:sizel
    t1{i} = rPla{i,1};
    % make a binary image where background signal is.
    thrs = graythresh(rPla{i,1})/3;
    tmask = im2bw(t1{i}, thrs * 1 );
    
    % calculate background level only inside the gonad using gonadal mask.
    % calculate background for individual z planes and save in 'mofm'
    tAll = t1{i}(tmask > 0);
    mofm(i) = mean(tAll)*1.5;     
end

% remove NaN from 'mofm' to calculate 'mmofm'
mofm(isnan(mofm)) = 1;
mom = mean(mofm);

mofm = mofm * mupr;


%%% threshold images and detect RNA spots
parfor i = 1:sizel 
% =============== peak detection method (1/3) ==========================
    [A, B] = FastPeakFind(t1{i}, mofm(i)*mupr);
    A = [A(1:2:end) A(2:2:end)];
        
    A(:,3) = 0;
    % remove peaks too close
    for j = 1:length(A(:,1))
        for k = j+1:length(A(:,1))
            dist = sqrt((A(k,1)-A(j,1))^2 + (A(k,2)-A(j,2))^2);
            if dist < pixr*2
                A(j,3) = 999;
                B(A(j,2),A(j,1)) = 0;
            end
        end
    end
    B = imdilate(B, se);
% ==================================================================

    temp = mean(t1{i}(B == 1));
    if isnan(temp)
        Pmofm(i) = 0;
    else
        Pmofm(i) = temp;
    end

%=========== visual detected RNAs =======
% figure,imshow(20*t1{i});
% figure,imshow(20*t1{i});
% hold on
% plot(A(:,1), A(:,2), 'r+');
% pause;
% close all
% end
% ---------------------------------------
    pMask{i} = B;
    
    fprintf('\nRNA%d: %d(th)/total %d images,... %d(th)/ %d z-planes.', g, f, lic, i, sizel);

end

Pmmofm = mean(Pmofm);
for i = 1:sizel
    fint{i} = t1{i} * Pmmofm  / Pmofm(i) * 2;
end


%%%%%%%%%%%%%%%%%%% 3D reconstitution %%%%%%%%%%%%%%%%%%%%%%%%%
%%% For mRNA (cytoplasmic spots)
tfrna = cat(3,pMask{:});
temp = bwconncomp(tfrna, 26);
conrna = regionprops(temp, 'Area', 'Centroid', 'BoundingBox', 'Image');
derna = struct2cell(conrna)';

% calculates intensity of detected blobs from original images.
iint = zeros(1,length(derna(:,1)));
for i=1:length(derna(:,1))
    iint(i) = 0;
    for j=1:derna{i,3}(6)
        cutimg = fint{round(derna{i,3}(3)),1};
        cutimg = cutimg(round(derna{i,3}(2)):round(derna{i,3}(2))...
            +derna{i,3}(5)-1, round(derna{i,3}(1)): round(derna{i,3}(1))+derna{i,3}(4)-1);

        iint(i) = iint(i) + sum(sum(derna{i,4}(:,:,j) .* double(cutimg)));
    end
end

derna = [derna num2cell(iint')];



if  isempty(derna) == 0
    if isempty(derna{1,1}) == 0
        % 'DETrnaP': | rna ID | x-size | y-size | z-size | 5th: total # pixel | 
        %           | total intensity | x centroid | y centroid | z
        %           centroid.
        DETrnaP = zeros(length(derna(:,1)),1);
        DETrnaP(:,1) = 1:length(derna(:,1));
        temp = cell2mat(derna(:,3));
        DETrnaP(:,2:4) = temp(:,4:6);
        temp = cell2mat(derna(:,2));
        DETrnaP(:,7:9) = temp(:,1:3);

        DETrnaP(:,6) = cell2mat(derna(:,5));
        DETrnaP(:,5) = cell2mat(derna(:,1));

        %%% remove false-positively detected blobs
        if isempty(DETrnaP) == 0
            % get mean of obj. size (pixels) and intensity (normalized: ind. int. - mean)
            mint = mean(DETrnaP(:,6));
            mpix = mean(DETrnaP(:,5));

            % take out objects detected less than 3 z-planes
            % take out objects dim/small.
            temp = DETrnaP(:,4) < 2 & DETrnaP(:,6) < mint*0.2 ;
            DETrnaP(temp,:) = [];
            
            
            if isempty(DETrnaP) == 0
                temp = DETrnaP(:,4) < 2 & DETrnaP(:,5) < mpix*0.2 ;
                DETrnaP(temp,:) = [];
            end
        end
    else
        DETrnaP = [];
    end
else
    DETrnaP = [];
end

% DETrna : |1-3: x,y,z-coordinates | 4: zero | 5: vol ratio to mean vol per spot |
%                   |6: total intensity| 7: z-plane of the center | 8: RNA
%                   ID
if ~isempty(DETrnaP)

    DETrna = DETrnaP(:,7:9);
    DETrna(:,3) = DETrna(:,3) * iminfo(7) / iminfo(6);
    DETrna(:,4) = 0;
    DETrna(:,5) = ceil(DETrnaP(:,5) / mean(DETrnaP(:,5))); %changed floor to ceil, deleted *2/3
%     DETrna(DETrna(:,5) == 0,5) = 1;
    DETrna(:,6) = DETrnaP(:,6);
    DETrna(:,7) = DETrnaP(:,9);
    DETrna(:,8) = DETrnaP(:,1);

    % remove mRNA with intensity 0
    DETrna(DETrna(:,6) <= 0,:) = [];
else
    DETrna = [];
end


