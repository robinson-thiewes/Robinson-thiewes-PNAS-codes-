%%% This code recognizes and records coordinates of nuclei and transcripts in each germline image in a folder.
function John_smFISH_v1_02_mpk1_v1_1( thresForExon, thresForIntron, thresForIntron2)
%% Clear workspace and create folders to save outputs.
clc; 
clearvars;
tic;

%%% Settings for image processing and analysis
%%% Create output folder and subfolders

%%%%%%------------- User Input (UI) system --------------------
% fpath = uigetdir('D:\smFISH images\');
% temp = regexp(fpath, '\', 'split');
% IMAGEfolderName = temp{end}; clear temp;
%-----------------------------------------------------------

%%%%%------------ Manual input -------------------------------
%% -------  for Windows
MASTER_Image_path = 'C:\Users\John\Documents\School\Kimble Lab\workable lif files on this PC'; % <<<<<-- Change the path!!!  Master image folder (one folder above the folder containing the image file to process).
IMAGEfolderName   =  'Test Data Sarah';           % <<<<<-- Change the path!!!  Input folder containing the image file (folder below the Master image folder).

MASTER_Output_path = 'C:\Users\John\Documents\School\Kimble Lab\Program Outputs'; % <<<<<-- Change the path!!!  Master Output folder
OUTPUTdir1=  'Test Data Sarah';
% OUTPUTdir1=  IMAGEfolderName;

%Create output folder and subfolders
MASTER_Image_path = fullfile(MASTER_Image_path, IMAGEfolderName, filesep);
cd(MASTER_Output_path); mkdir(OUTPUTdir1); cd(OUTPUTdir1);
savepath = fullfile(MASTER_Output_path,OUTPUTdir1, filesep);
%%-----------------------------------------------------------

% % %%% --------  for MAC --------------------------------------
% MASTER_Image_path = strcat('/Users/sarahrobinson/Desktop/smFISH images folder/'); % <<<<<-- Change the path!!!  Master image folder (one folder above the folder containing the image file to process).
% IMAGEfolderName   =  'N2';                   % <<<<<-- Change the path!!!  Input folder containing the image file (folder below the Master image folder).
%
% MASTER_Output_path=  '/Users/sarahrobinson/Desktop/smFISH images folder/';              % <<<<<-- Change the path!!!  Output folder
% %   OUTPUTdir1=  'N2';
%
% %  Create output folder and subfolders
% MASTER_Image_path = strcat(MASTER_Image_path, IMAGEfolderName, '/');
% cd(MASTER_Output_path); mkdir(OUTPUTdir1); cd(OUTPUTdir1);
% savepath = strcat(MASTER_Output_path,OUTPUTdir1,'/');
% % %%%-----------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Threshold for image processing (-0.5 to inf; normally 0 - 1). Change parameter if not detect all spots
% Thrshold for Exon/Intron:
% N2: .5/.5     lag-2(q411)/+: .5/.7     glp-1(q46)/+: .5/.7     lag-1(q385)/+: .3/.7   JK5008: .65/.45 .
% lag-3(ok387) het: .5/.4
% q224:
% Erika's mex-3: .5/.5
thresForExon =         0.65     ;
thresForIntron =       0.6    ;
thresForNuc =    .5      ;
%nthfile =         1     ;

% Designate the channel for RNA.
% Put   '1' for intron,     '2' for exon,      '3' for else.
ChOrder = [     1       2        ];


radius =     2.5     ;   % set the Radius of ROI as desired (um)
vLim =      3        ;   % maxium distance for Voronoi cell.
nrange = [   1.2     2.8   ];    % define range of nuclear radius in um.
sensi =     0.96   ;     % sensitivity for nuclear circle detection (higher: more flex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\n\t\tThe radius of ROI is %2.1f um.\n\n', radius);

[nbf, listSort] = read_lif_files(MASTER_Image_path);

af = cell(999,10,nbf);

%% Process multiple lif files sequentially
for h=1:nbf    % h is a counter for lif files (with diff conditions) in the folder.
    %     clearvars -except af pix sensi nrange st_names indeX liflist h listSort pDir1...
    %         stname_forNum CIpercent nbf radius label1 label2 label3 label4 fpath savepath...
    %         thresForExon thresForIntron thresForNuc ChOrder vLim kkk kno nthfile ifo;
    
    %%% import lif files in the same conditions (e.g. strain, temp)
    lifdat = cell(1,4);
    fToRead = strcat(MASTER_Image_path, listSort{h,1});
    
    % read in total # of image stacks
    lifRead = bfGetReader(fToRead);
    omeMeta = lifRead.getMetadataStore();
    lic = omeMeta.getImageCount();
    
    
    %% Process and analyze individual images one by one.
    %     kk = 4:5;
    
    for f =  1:lic     % open f-th image-stack in the lif file.
    % %
    %         if f == 1
    %             continue
    %         end
    % %
    if f == 3
        thresForNuc =    .8      ;
    end
    %         elseif f == 1
    %             thresForIntron =    0.3      ;
    %
    %         end
    
    lifdat = bfopenOne(fToRead,f);
    
    %%% read metadata
    meta = lifdat{1, 4};
    Nch = meta.getPixelsSizeC(f-1).getValue();
    voxelSizeXdefaultValue = meta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
    voxelSizeXdefaultUnit = meta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
    voxelSizeX = meta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
    voxelSizeY = meta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER);
    voxelSizeZ = meta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER);
    
    iminfo = [
        Nch % # channel
        meta.getPixelsSizeX(f-1).getValue(); % image width, pixels
        meta.getPixelsSizeY(f-1).getValue(); % image height, pixels
        meta.getPixelsSizeZ(f-1).getValue(); % number of Z slices
        voxelSizeX.doubleValue();  % in µm
        voxelSizeY.doubleValue();  % in µm
        voxelSizeZ.doubleValue()   % in µm
        ];
    
    
    %%% in case above method does not work for reading # of z-planes-----
    %         temp = meta.getPixelsPhysicalSizeZ(f-1).getValue(); % in um, works on ver.R2014b and +.
    %         if isnan(temp) == 1
    %             str2double(lifdat{1,2}.get('ATLConfocalSettingDefinition|Quantity|Value 0')) * 1000000;
    %         end
    %         iminfo = [iminfo; temp];
    %---------------------------------------------------------------------
    
    %%% Import image information from metadata
    % find Nuc, phase or trx channel by color (channel label does not work).
    % Chinfo: |col 1: color info in the order of acquisition|
    Chinfo = zeros(Nch,2);
    for i = 1:Nch
        Chinfo(i,1) = meta.getChannelColor(f-1,i-1).getValue();
    end
    
    % Nuc | phase | Yellow (often 561/594) | Red (often 633/594) | green (often LMN-1 or other Ab with Alexa488)
    % 1: DAPI(blue or cyan), 2: Phase, 3: Yellow, 4: Red, 5: Magenta, 6: Green
    Chinfo(Chinfo(:,1) == 16777215 | Chinfo(:,1) == 65535,2) = 1; % DAPI, blue or cyan
    Chinfo(Chinfo(:,1) == -1,2) = 2; % gray, NOT phase channel
    Chinfo(Chinfo(:,1) == -65281,2) = 3; % yellow
    Chinfo(Chinfo(:,1) == -16711681,2) = 4; % magenta
    Chinfo(Chinfo(:,1) == -16776961,2) = 5; % red
    Chinfo(Chinfo(:,1) == 16711935,2) = 6; % green
    if find (Chinfo(:,1) == -1, 1, 'last') == Nch
        Chinfo(Nch,2) = 9; % phase
    end
    
    loc = cell2mat(strfind(lifdat{1,1}(:,2), 'C='))+2;
    chameta = lifdat{1,1}(:,2);
    leng = length(chameta(:,1));
    chainfo = zeros(leng,1);
    parfor i = 1:leng
        temp = char(chameta(i,1));
        chainfo(i) = str2double(temp(loc(i)));
    end
    
    
    
    
    
    %% Detect a germline outline using ATS channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Seperate planes
    % The last number on line below determines which color channel to use.
    PhaPla = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == 6), :);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %changed to 4 from 6 (green) above
    pix = round( 0.22 / mean(iminfo(5:6,1)));
    if pix == 0
        pix = 1;
    end
    
    zproj = cat(3, PhaPla{:,1});
    mip = max(zproj, [], 3);
    mipBW = imbinarize(mip, graythresh(mip)/3);
    
    
    %%% "cyl_axis" contains: (unit: pixels)
    %%% (1) x-coor | (2) y-coor | (3) z-coor | (4) radius of the cell on that point.|
    [cyl_axis, germlineBound,cmask] = DefOutlineATS(mip, pix, iminfo);
    
    %%% align germline orientation (distal end to the left)
    if mean(cyl_axis(1:10,4)) > mean(cyl_axis(end-9:end,4))
        cyl_axis(:,1) = iminfo(2)-cyl_axis(:,1)+1;
        cmask = flip(cmask, 2);
        %lifdat is {1x4}, {1,1} is image pixels: {1,1}{:,1}, image name: {1,1}{:,2}
        lifdat{1,1}(:,1) = cellfun(@(x) flip(x,2), lifdat{1,1}(:,1), 'uni', 0); %flip images if germline is flipped
    end
    
    
    %% Detect Nuclei
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % seperate Nuc channel to detect nuclei.
    NucPla = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == 1), :);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Detect nuclei and remove one nucleus if two are overlapped.
    % NDNuc saves not overlapped Nuclei that are in DetNuc (detected Nuc).
    tnran = round(nrange ./ iminfo(5:6,1)');
    tol = round( 1.2 / mean(iminfo(5:6,1))); %%%%% tolerance for overlapped circles (um).
    
    % Detect nuclei from image slices and make nuclear sphere.
    
    
    [Nuc, thrNuc, aftNuc, tickn] = DetectNucleus(NucPla, iminfo, pix, tnran, sensi, tol, f, lic, thresForNuc);
    
    % Detect nuclei using bwconncomp and cross-validate the result with nuclei detected by hough transf.
    % 'nucs' : | x-coor | y | z | radius | ave. circularity |
    %           | std. circularity | total DAPI sig. area (# pixel) |
    nucs = DetectNucBlobs(thrNuc, aftNuc, Nuc, iminfo, tnran);
    
    % Remove nuclei that are outside the gonad
    nucs(:,8) = 1;
    for i = 1:size(nucs,1)
        r1 = round(nucs(i,2))-3;
        r2 = round(nucs(i,2))+3;
        r3 = round(nucs(i,1))-3;
        r4 = round(nucs(i,1))+3;
        if r1 < 1
            r1 = 1;
        end
        if r2 > iminfo(3)
            r2 = iminfo(3);
        end
        if r3 < 1
            r3 = 1;
        end
        if r4 > iminfo(2)
            r4 = iminfo(2);
        end
        if sum(sum(mipBW(r1:r2, r3:r4))) == 0
            nucs(i,8) = 0;
        end
    end
    nucs(nucs(:,8) == 0,:) = [];
    nucs(:,8) = [];
    
    % Visualize of nuclei
    %----------- max projection  ------------------
    zproj = cat(3, NucPla{:,1});
    mipNuc = max(zproj, [], 3);
    %         figure, imshow(mip);
    
    
    %% Detect RNA spots
    % Process different RNA (channels) separately.
    rch = [3 4 5 6];
    chname = {'yel', 'mag', 'red', 'grn'};
    sele = ismember(rch, Chinfo(:,2));
    chnn = rch(sele); %chnn contains the RNA channels
    chn = Chinfo(ismember(chnn,Chinfo),2)';
    nch = length(chn);
    chname = chname(sele);
    prna = cell(nch,1);
    nuc = cell(nch,1);
    
    for i = 1:nch
        rnaPla.(chname{i}) = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == chn(i)), :);
        if i == 1
            k = 10; % 1st RNA (channel) image is saved in 10th column
        elseif i == 2
            k = 11; % 2nd RNA (channel) image is saved in 11th column
        else
            k = 12; % 3rd RNA (channel) image is saved in 12th column
        end
        zproj = cat(3, rnaPla.(chname{i}){:,1});
        zproj = max(zproj, [], 3);
        af{f,k} = zproj;
    end
    
    %%% process channels seperately for detecting RNA
    % 'blobrna': | rna ID | x-size | y-size | z-size | 5th: total # pixel |
    %            | total intensity | x centroid | y centroid | z centroid.
    
    % blrna: |Area | Centroid | BoundingBox | Image (stack) | total fluor. intensity.
    %        | 6th: matching nuc #(row # in 'nuc')
    
    % 'rnas' 1:3 col = x,y,z coor, 4 col = nuc assciated? | 5th: vol ratio to mean vol per spot |
    %               | 6th: total sig intensity (trx activity) | 7: z-coor in plane # | 8: RNA ID.
    %               | 9: count for Voronoi matching | 10: Voronoi count w/ limit (vLim) |
    %               | 11: # total ROI (overlapped) | 12: distance to closest nucleus|
    
    % 'nuc for ATS' : | x-coor | y | z | radius | ave. circularity |
    %         | 6th: std. circularity | total DAPI sig.  of all nuclei | # trx sites | RNA ID .
    
    % 'nuc for mRNA' : | 1:3 xyz-coor | 4: radius | 5: ave. circularity |
    %           | 6: std. circularity | 7: total DAPI sig. area (# pixel) | 8: nuc ID.
    %           | 9: Voronoi mRNA count | 10: Voronoi count w/ limit (vLim, 3um) | 11: # mRNA in ROI w/ radius (2.5um)  |
    %  (aft_analysis->)       | 12th col: dist betw. nuc & cell outline | 13: 0 = cells on cortex, 1 = cells in rachis |
    
    % 'ExATSnuc' | 1st: x-coor | y | z | 4th: radius | 5th: mean circularity |
    %            | 6th: std circularity nuc | 7th: total intensity of DAPI |
    %            | 8th: # of ATS per nuc | 9th: summed intensity of dectected ATS |
    
    
    % 'ExATSrna'  | 1-3: x, y, z coor | 4: nuc ID | 5th: estimated true ATS within spot |
    %             | 6th: summed intensity of detection | 7th: z plane | 8th: mRNA ID |
    %
    
    pixr = round( 0.1 / mean(iminfo(5:6,1)));
    
    parfor g = 1:nch
        rPla = rnaPla.(chname{g})(:,1);
        %%%%%%%%%%%%%% Detect RNA spots %%%%%%%%%%%%%%%%%%%%%%%%%
        %-------------- process INTRON -----------------------------
        if ChOrder(g) == 1    % 1 is for processing intron images.
            % mom1 or mom2 is the mean background signal intensity from all z-planes in a gonad.
            [blobrna, blrna, mom1(g)] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, lic, g, thresForIntron);
            
            if isempty(blobrna) == 0
                % determine which nuclei RNA spots belong.
                [prnas, nuct] = MatchrnaNuc(thrNuc, blobrna, nucs, iminfo, pixr, 1.2);
                
                % Put results of RNA1 and RNA2 in separate struct.
                prna{g} = prnas;
                nuc{g} = nuct;
            else
                prna{g} = [];
                nuc{g} = [];
            end
            
            % ------------- process EXON image ----------------------------
        elseif ChOrder(g) == 2   % 2 is for processing exon images.
            [DETrna, derna, mom2(g)] = DetectRNAExon(rPla, pixr, iminfo, f, lic, g, thresForExon);
            
            if isempty(DETrna) == 0
                
                % Put results of RNA1 and RNA2 in separate struct.
                DETrna(:,8) = 1:length(DETrna(:,1)); % mRNA ID.
                [mrnas, rnas_atsEx, nucm, nuc_atsEx] = MatchmrnaNuc(DETrna, thrNuc, nucs, iminfo, radius, vLim, pixr);
                
                prna{g} = mrnas; % save info for detected mRNA
                prna{g}(:,4) = 1;
                nuc{g} = nucm;
                
                %%% record detected ATS from exon channel
                ExATSrna{g} = rnas_atsEx;
                ExATSnuc{g} = nuc_atsEx;
                
            else
                prna{g} = [];
                nuc{g} = [];
            end
        elseif ChOrder(chn(g)-2) == 3   % 3 is for processing etc.
            
        end
        
    end
    
    
    
    %% cross-validate trx sites by Intron probe and by Exon probe.
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
                if i == 14
                    i = 14;
                end
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
                        ExATSrna_ATS(i,9) = 2121
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
    
    
    
    %% Make summary data
    % 'af': | 1st: nuc1 info | 2nd: rna1 info (yellow) | 3rd: nuc2 info | 4th: rna2 info (red) | .
    %           | 5th: cyl_axis info | 6th: cell outline mask | 7th: z-prj nuclear channel |.
    %           | 8th: mean bg intensity from 1st channel | 9th: mean bg. noise from 2nd channel|
    %           | 10th: z-projected RNA image (1st channel) | 11th: z-prj. 2nd RNA channel |
    %  (aft_analysis->)          | 12th col: distance distal tip & first chain nuc | 12: # nuc in rachis |
    for i = 1:nch
        if ~isempty(prna{i})
            temp = prna{i}(prna{i}(:,4) ~= 0,:);
            af(f,i*2-1, nbf) = nuc(i);
            af(f,i*2, nbf) = {temp};
        else
            af(f,i*2-1, nbf) = {[]};
            af(f,i*2, nbf) = {[]};
        end
    end
    
    af{f,5, nbf} = cyl_axis;
    af{f,6, nbf} = cmask;
    af{f,7, nbf} = mip;
    af{f,12,nbf} = ExATSrna_ATS;
    if exist('mom1', 'var')
        af{f,8, nbf} = mom1(mom1~=0);
    end
    if exist('mom2', 'var')
        af{f,9, nbf} = mom2(mom2~=0);
    end
    
    %%% record information for nuc (12 col) & detected ATS (13th col) from exon channel in 'af'
    if exist('ExATSrna', 'var')
        for i = 1:size(ExATSrna,2) %changed 1 to 2 because exon 2nd channel
            if ~isempty(ExATSrna{i})
                temp = ExATSrna{i}(ExATSrna{i}(:,4) ~= 0,:);
                af(f,i*2-1 +11, nbf) = ExATSnuc(i);
                af(f,i*2 +11, nbf) = {temp};
            else
                af(f,i*2-1 +11, nbf) = {[]};
                af(f,i*2 +11, nbf) = {[]};
            end
        end
    end
    
    cd(savepath);
    save(strcat('smFISHworkspace_middle_',num2str(f)), '-regexp', '^(?!(lifdat|lifRead|meta|omeMeta)$).');
    
    end
end

if nbf == 1
    af(lic+1:end,:) = [];
else
    af(lic+1:end,:,nbf) = [];
end

cd(savepath);
save('smFISHworkspace_final', '-regexp', '^(?!(lifdat|lifRead|meta|omeMeta)$).');

eTime = toc;
fprintf('\n\t\tElapsed time is %2.1f minutes.\n\n', eTime/60);
end




function [nbf, liflist] = read_lif_files(MASTER_Image_path)
    %% Read lif files
    %%% lif file format: separate diff condition by ' ' (space). Order or # of
    % statements does not matter. e.g. 040914 rrf-1 emptyRNAi_na lag-3_na.lif.
    % liflist = row1: lif file name w/o '.lif'; row2: name with sorted; conditions; row3: strain numbering.
    [liflist] = dir(strcat(MASTER_Image_path,'*.lif'));
    liflist = liflist(~[liflist.isdir]);
    liflist = {liflist.name}';
    nbf = numel(liflist);

end

