%%% This code recognizes and records coordinates of nuclei and transcripts in each germline image in a folder. 
% splitProbe_v3_01('Split probe mpk1 output',.5,.5,.5)
% v3: Original Version
% v3_01: Copy of original version
% v3_02: validate as separate function
   function splitProbe_v3_01_2(OUTPUTdir1, thresForExon, thresForIntron1, thresForIntron2)
%% Clear workspace and create folders to save outputs.
% clc; clearvars; 
tic; %times how long it takes to run program

%%% Settings for image processing and analysis
%%% Create output folder and subfolders

%%%%%%------------- User Input (UI) system --------------------
% fpath = uigetdir('D:\smFISH images\');
% temp = regexp(fpath, '\', 'split');
% IMAGEfolderName = temp{end}; clear temp;
%-----------------------------------------------------------

%%%%%%------------ Manual input -------------------------------
%%% -------  for Windows
MASTER_Image_path = 'C:\Users\John\Documents\School\Kimble Lab\workable lif files on this PC\'; % <<<<<-- Change the path!!!  Master image folder (one folder above the folder containing the image file to process).
IMAGEfolderName   =  'mpk1 split probe N2';           % <<<<<-- Change the path!!!  Input folder containing the image file (folder below the Master image folder).

MASTER_Output_path = 'C:\Users\John\Documents\School\Kimble Lab\Program Outputs\'; % <<<<<-- Change the path!!!  Master Output folder
OUTPUTdir1=  'Split probe mpk1 output';
% OUTPUTdir1=  IMAGEfolderName;

%Create output folder and subfolders
MASTER_Image_path = strcat(MASTER_Image_path, IMAGEfolderName, '\');
cd(MASTER_Output_path); mkdir(OUTPUTdir1); cd(OUTPUTdir1);
savepath = strcat(MASTER_Output_path,OUTPUTdir1,'\');
%%%-----------------------------------------------------------
                                                     
% %%% --------  for MAC --------------------------------------                                                                                         
% MASTER_Image_path = strcat('/Users/sarahrobinson/Desktop/smFISH images folder/'); % <<<<<-- Change the path!!!  Master image folder (one folder above the folder containing the image file to process).
%  IMAGEfolderName   =  'N2';                   % <<<<<-- Change the path!!!  Input folder containing the image file (folder below the Master image folder).
%  
%  MASTER_Output_path=  '/Users/sarahrobinson/Desktop/smFISH images folder/';              % <<<<<-- Change the path!!!  Output folder
% %  OUTPUTdir1=  'N2';
%  
% %  Create output folder and subfolders
%  MASTER_Image_path = strcat(MASTER_Image_path, IMAGEfolderName, '/');
%  cd(MASTER_Output_path); mkdir(OUTPUTdir1); cd(OUTPUTdir1);
%  savepath = strcat(MASTER_Output_path,OUTPUTdir1,'/');
% %%%-----------------------------------------------------------  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Threshold for image processing (-0.5 to inf; normally 0 - 1). Change parameter if not detect all spots
% Thrshold for Exon/Intron:
% N2: .5/.5     lag-2(q411)/+: .5/.7     glp-1(q46)/+: .5/.7     lag-1(q385)/+: .3/.7   JK5008: .65/.45 .
% lag-3(ok387) het: .5/.4
% q224: 
% Erika's mex-3: .5/.5
% thresForExon =         0.5     ; 
% thresForIntron1 =       0.5     ;

thresForNuc =    .5      ; %threshold for nuclei
nthfile =         1     ;    %labels file number   

% Designate the channel for RNA.
% Put   '1' for intron first half,     '2' for exon,      '3' for intron second half probe.
ChOrder = [    3    1    2     ];


radius =     2.5     ;   % set the Radius of ROI as desired (um)
vLim =      3        ;   % maxium distance for Voronoi cell.
nrange = [   1.2     2.8   ];    % define range of nuclear radius in um.
sensi =     0.96   ;     % sensitivity for nuclear circle detection (higher: more flex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%% Open matlabpool for parallel computing
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end

fprintf('\n\t\tThe radius of ROI is %2.1f um.\n\n', radius); 


%% Read lif files 
%%% lif file format: separate diff condition by ' ' (space). Order or # of
% statements does not matter. e.g. 040914 rrf-1 emptyRNAi_na lag-3_na.lif.
% liflist = row1: lif file name w/o '.lif'; row2: name with sorted; conditions; row3: strain numbering.
[liflist] = dir(strcat(MASTER_Image_path,'*.lif')); %opens directory for images
liflist = liflist(~[liflist.isdir]);
liflist = {liflist.name}';
nbf = numel(liflist); %makes array element based on liflist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Choose n-th lif file to process
nbf = 1;
liflisttemp = liflist(           nthfile           );
liflist=liflisttemp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% gather lif files with the same labels (experimental conditions).
liflist = strrep(liflist,'.lif','');
lifSplit = regexp(liflist, ' ', 'split');

parfor i = 1:nbf %set up parallel processing
    lifSplit{i} = lifSplit{i}(cellfun('isempty', regexp(lifSplit{i},'^-?\d+$')));
    lifSplit{i} = sort(lifSplit{i}); 
    liflist{i,2} = strjoin(lifSplit{i},'_');
end

[listSort,indeX] = sortrows(liflist, 2);
listSort{1,3} = 1; liflist{indeX(1),3} = 1;

count = 1;
st_names = cell(nbf,1);
st_names{1} = listSort{1,2};
for i = 1:nbf-1
    if strcmp(listSort{i+1,2},listSort{i,2}) == 1
        listSort{i+1,3} = listSort{i,3};
        liflist{indeX(i+1),3} = liflist{indeX(i),3};
    else
        count = count+1;
        listSort{i+1,3} = count;
        liflist{indeX(i+1),3} = count;
        st_names{count,1} = listSort{i+1,2};
    end
end

if nbf ~= 1
    st_names(count+1,:) = [];
end

clear count
st = length(st_names);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stname_forNum = cell(1,listSort{end,end}); 
af = cell(999,17,nbf); %premakes af cell to be filled in later


%% Process multiple lif files sequentially
for h=1:nbf    % h is a counter for lif files (with diff conditions) in the folder.
%     clearvars -except af pix sensi nrange st_names indeX liflist h listSort pDir1...
%         stname_forNum CIpercent nbf radius label1 label2 label3 label4 fpath savepath...
%         thresForExon thresForIntron thresForNuc ChOrder vLim kkk kno nthfile ifo;
        
%%% import lif files in the same conditions (e.g. strain, temp)
lifdat = cell(1,4);
fToRead = strcat(MASTER_Image_path, listSort{h,1}, '.lif');

% read in total # of image stacks
lifRead = bfGetReader(fToRead);
omeMeta = lifRead.getMetadataStore();
lic = omeMeta.getImageCount();


    %% Process and analyze individual images one by one.
%     kk = 4:5;
    for f =  1:lic     % open f-th image-stack in the lif file.
% %         
        if f == 14 || f == 15 || f== 16 || f== 17 || f == 2 || f == 9 %(||) means or
            continue %2 and 9 cause problems
        end %this is for skipping images, f==image number
%         
%         if f == 7
%             thresForNuc =    .8      ;
%         elseif f == 5
%             thresForNuc =    1      ;
%         else
%             thresForNuc =    .5      ;
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
        Chinfo(Chinfo(:,1) == -65281,2) = 3; % yellow, 1st RNA
        Chinfo(Chinfo(:,1) == -16711681,2) = 4; % magenta, 2nd RNA?
        Chinfo(Chinfo(:,1) == -16776961,2) = 5; % red, 3rd RNA     
        Chinfo(Chinfo(:,1) == 16711935,2) = 6; % green, Ab staining? (e.g. Alexa-488)
        if find (Chinfo(:,1) == -1, 1, 'last') == Nch
            Chinfo(Nch,2) = 9; % phase
        end

        loc = cell2mat(strfind(lifdat{1,1}(:,2), 'C='))+2;
        chameta = lifdat{1,1}(:,2);
        leng = length(chameta(:,1));
        chainfo = zeros(leng,1);
        parfor i = 1:leng %parallel processing
            temp = char(chameta(i,1));
            chainfo(i) = str2double(temp(loc(i)));
        end



        
 %% Detect a germline outline using ATS channel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Seperate planes
        % The last number on line below determines which color channel to use.
        PhaPla = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == 6), :);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%changed from 4 (magenta) to 6 (green) above
        pix = round( 0.22 / mean(iminfo(5:6,1))); 
        if pix == 0
            pix = 1;
        end

        zproj = cat(3, PhaPla{:,1});
        mip = max(zproj, [], 3);
        mipBW = imbinarize(mip, graythresh(mip)/3); %changed im2bw to imbinarize


        %%% "cyl_axis" contains: (unit: pixels)
        %%% (1) x-coor | (2) y-coor | (3) z-coor | (4) radius of the cell on that point.| 
        [cyl_axis, germlineBound,cmask] = DefOutlineATS(mip, pix, iminfo);

        %%% align germline orientation (distal end to the left)
        if mean(cyl_axis(1:10,4)) > mean(cyl_axis(end-9:end,4))
            cyl_axis(:,1) = iminfo(2)-cyl_axis(:,1)+1;
            cmask = flip(cmask, 2);
            lifdat{1,1}(:,1) = cellfun(@(x) flip(x,2), lifdat{1,1}(:,1), 'uni', 0);
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
%         zproj = cat(3, NucPla{:,1});
%         mip = max(zproj, [], 3);
%         figure, imshow(mip);


        %% Detect RNA spots
        % Process different RNA (channels) separately.
        rch = [3 4 5 6];
        chname = {'yel', 'mag', 'red', 'grn'};
        sele = ismember(rch, Chinfo(:,2));
        chnn = rch(sele); 
        chn = Chinfo(ismember(chnn,Chinfo),2)';
        nch = length(chn);
        chname = chname(sele);
        prna = cell(nch,1);
        nuc = cell(nch,1);

        for i = 1:nch
            rnaPla.(chname{i}) = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == chn(i)), :);
            if i == 1
                k = 13; % 1st RNA (channel) image is saved in 13th column
            elseif i == 2
                k = 14; % 2nd RNA (channel) image is saved in 14th column
            elseif i == 3 %changed else ... line
                k = 15; % 3rd RNA (channel) image is saved in 15th column
            end
            
            
            zproj = cat(3, rnaPla.(chname{i}){:,1});
            zproj = max(zproj, [], 3);
            af{f,k} = zproj;
        end
%%  %%%%%%%%%%%%%% Detect RNA spots %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% process channels seperately for detecting RNA
        % 'blobrna': | 1st: rna ID | 2nd: x-size 
        %            | 3rd: y-size | 4th: z-size 
        %            | 5th: total # pixel | 6th: total intensity
        %            | 7th: x centroid | 8th: y centroid 
        %            | 9th: z centroid.

        % blrna: | 1st: Area | 2nd: Centroid 
        %        | 3rd: BoundingBox | 4th: Image (stack) 
        %        | 5th: total fluor. intensity.
        %        | 6th: matching nuc #(row # in 'nuc') 

        % 'rnas' | 1st: x coordinate | 2nd: y coordinate
        %        | 3rd: z coordinate | 4th: nuc assciated? 
        %        | 5th: vol ratio to mean vol per spot 
        %        | 6th: total sig intensity (trx activity) 
        %        | 7th: z-coor in plane # | 8th: RNA ID.
        %        | 9th: count for Voronoi matching | 10th: Voronoi count w/ limit (vLim) 
        %        | 11th: # total ROI (overlapped) | 12th: distance to closest nucleus

        % 'nuc for ATS': | 1st: x coordinate | 2nd: y coordinate
        %                | 3rd: z coordinate | 4th: radius 
        %                | 5th: ave. circularity | 6th: std. circularity
        %                | 7th: total DAPI sig.  of all nuclei | 8th: # trx sites 
        %                | 9th: RNA ID .

        % 'nuc for mRNA': | 1st: x coordinate | 2nd: y coordinate
        %                 | 3rd: z coordinate | 4th: radius 
        %                 | 5th: ave. circularity | 6th: std. circularity
        %                 | 7th: total DAPI sig. area (# pixel) | 8th: nuc ID.
        %                 | 9th: Voronoi mRNA count | 10th: Voronoi count w/ limit (vLim, 3um) 
        %                 | 11th: # mRNA in ROI w/ radius (2.5um) (aft_analysis->)
        %                 | 12th col: dist betw. nuc & cell outline | 13th: 0 = cells on cortex, 1 = cells in rachis 

        pixr = round( 0.1 / mean(iminfo(5:6,1)));

        parfor g = 1:nch
            rPla = rnaPla.(chname{g})(:,1);   
           
            %-------------- process INTRON first half -----------------------------
            %This g mean starting with position 1 to total number of total channels (3 in this case).
            %Looks to see which number in channel order is logical =1
            if ChOrder(g) == 1    % 1 is for processing intron 1st half images.
                % mom1, mom2, mom3 is the mean background signal intensity from all z-planes in a gonad.
                [blobrna, blrna, mom1(g)] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, lic, g, thresForIntron1);

                blobrna1 = blobrna
                blrna1 = blrna
                
                if isempty(blobrna1) == 0
                    % determine which nuclei RNA spots belong. 
                    [prnas, nuct] = MatchrnaNuc_2(thrNuc, blobrna1, nucs, iminfo, pixr, 1.2);

                    prnas1 = prnas
                    nuct1 = nuct
                    
                    % Put results of RNA1, RNA2, and RNA3 in separate struct.
                    prna{g} = prnas1;
                    nuc{g} = nuct1;
                else
                    prna{g} = [];
                    nuc{g} = [];
                end
                
            
            %-------------- process INTRON second half -----------------------------
            elseif ChOrder(g) == 3    % 3 is for processing intron second half images.
                % mom3 is the mean background signal intensity from all z-planes in a gonad.
                [blobrna, blrna, mom3(g)] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, lic, g, thresForIntron2);

                blobrna3=blobrna
                blrna3=blrna
                
                if isempty(blobrna3) == 0
                    % determine which nuclei RNA spots belong. 
                    [prnas, nuct] = MatchrnaNuc(thrNuc, blobrna3, nucs, iminfo, pixr, 1.2);

                    prnas3 = prnas
                    nuct3 = nuct
                    
                    % Put results of RNA1 and RNA2 in separate struct.
                    prna{g} = prnas3;
                    nuc{g} = nuct3;
                else
                    prna{g} = [];
                    nuc{g} = [];
                end

            % ------------- process EXON image ----------------------------    
            elseif ChOrder(g) == 2   % 2 is for processing exon images.
                
                %added to see if exon can be used for ATS (4 used to designate exon ATS)
%                 [blobrna, blrna, mom4(g)] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, lic, g, thresForExon);
% 
%                 blobrna4=blobrna
%                 blrna4=blrna
                
%                 if isempty(blobrna4) == 0
%                     % determine which nuclei RNA spots belong. 
%                     [prnas, nuct] = MatchrnaNuc(thrNuc, blobrna4, nucs, iminfo, pixr, 1.2);
% 
%                     prnas4 = prnas
%                     nuct4 = nuct
                    
%                     % Put results of exon ATS and nuclear detection in separate arrays.  Needed to use g+2 so that data not overwritten

%                     prna{g+2} = prnas4;
%                     nuc{g+2} = nuct4;
%                 else
%                     prna{g+2} = [];
%                     nuc{g+2} = [];
%                 end
                
                %RNA detection to include mRNAs
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
%             elseif ChOrder(chn(g)-2) == 3   % 3 is for processing etc.
% 
            end

        end 

        % 'prna' | 1st: x-coor  | 2nd: y-coor
        %        | 3rd: z-coor  | 4th: nuc #
        %        | 5th: 1 or 2 spots   | 6th: total signal intensity
        %        | 7th: z-coor in plane #
        
        % 'mrnas' | 1st: x-coor  | 2nd: y-coor
        %         | 3rd: z-coor  | 4th: nuc associated or cell associated
        %         | 5th: volume ratio to mean vol/spot   | 6th: total signal intensity
        %         | 7th: z-coor in plane #   | 8th: RNA ID
        %         | 9th: count for Voronoi match    | 10th: Voronoi count with limit
        %         | 11th: # total ROI        | 12th: distance to closest nuc
        %         
 
        %% cross-validate trx sites by Intron probe and by Exon probe.
    [prna, nuc] = splitprobe_validate_v3_02(prna, nuc, ChOrder, iminfo);

        %% Make summary data
        % 'af': | 1st: nuc1 info | 2nd: rna1 info  
%               | 3rd: nuc2 info | 4th: rna2 info 
%               | 4th: nuc3 info | 5th: rna3 info 
%               | 7th: cyl_axis info | 8th: cell outline mask 
%               | 9th: z-prj nuclear channel | 10th: mean bg intensity from 1st channel 
%               | 11th: mean bg. noise from 2nd channel| 12th: mean bg. intensity from 3rd channel 
%               | 13th: z-projected RNA image (1st channel)| 14th:z-prj. 2nd RNA channel
%               | 15th: z-projected RNA image (3rd channel)| 16th: nuc exon detected
%               | 17th: exon dectected ATS
        
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

        af{f,7, nbf} = cyl_axis;
        af{f,8, nbf} = cmask;
        af{f,9, nbf} = mip;
        if exist('mom1', 'var')
            af{f,10, nbf} = mom1(mom1~=0);
        end
        if exist('mom2', 'var')
            af{f,11, nbf} = mom2(mom2~=0);
        end
        if exist('mom3', 'var')
            af{f,12, nbf} = mom3(mom3~=0);
        end

        %%% record information for nuc (12 col) & detected ATS (13th col) from exon channel in 'af'
        for i = 1:size(ExATSrna,1)
            if ~isempty(ExATSrna{i})
                temp = ExATSrna{i}(ExATSrna{i}(:,4) ~= 0,:);
                af(f,i*2-1 +15, nbf) = ExATSnuc(i);
                af(f,i*2 +15, nbf) = {temp};
            else
                af(f,i*2-1 +15, nbf) = {[]};
                af(f,i*2 +15, nbf) = {[]};
            end
        end
        
        cd(savepath);
        save(strcat('smFISHworkspace_middle_',num2str(f)), '-regexp', '^(?!(lifdat|lifRead|meta|omeMeta)$).');
    end

    if nbf == 1
        af(lic+1:end,:) = [];
    else
        af(lic+1:end,:,nbf) = [];
    end
       
    
end


cd(savepath); %saves file in designated spot in beginning
save('smFISHworkspace_final', '-regexp', '^(?!(lifdat|lifRead|meta|omeMeta)$).');

eTime = toc;
fprintf('\n\t\tElapsed time is %2.1f minutes.\n\n', eTime/60); 