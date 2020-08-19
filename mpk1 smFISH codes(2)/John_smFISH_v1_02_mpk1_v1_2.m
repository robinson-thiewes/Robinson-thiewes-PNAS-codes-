%%% This code recognizes and records coordinates of nuclei and transcripts in each germline image in a folder.
function John_smFISH_v1_02_mpk1_v1_2() %thresForExon, thresForIntron, thresForIntron2)
%% Clear workspace and create folders to save outputs.
clc; clearvars; tic;

%%%%%------------ Manual input -------------------------------
MASTER_Image_path = 'D:\smFISH images'; % <<<<<-- Change the path!!!  Master image folder (one folder above the folder containing the image file to process).
IMAGEfolderName   =  'images';           % <<<<<-- Change the path!!!  Input folder containing the image file (folder below the Master image folder).

MASTER_Output_path = 'D:\smFISH images'; % <<<<<-- Change the path!!!  Master Output folder
OUTPUTdir1=  'Matlab processing';

MASTER_Image_path = fullfile(MASTER_Image_path, IMAGEfolderName, filesep);
cd(MASTER_Output_path); mkdir(OUTPUTdir1); cd(OUTPUTdir1);
savepath = fullfile(MASTER_Output_path,OUTPUTdir1, filesep);
%%%%%-----------------------------------------------------------



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
outline_channel = 6 ;
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
    %%% import lif files in the same conditions (e.g. strain, temp)
    fToRead = strcat(MASTER_Image_path, listSort{h,1});
    % read in total # of image stacks
    lifRead = bfGetReader(fToRead);
    omeMeta = lifRead.getMetadataStore();
    lic = omeMeta.getImageCount();
    %% Process and analyze individual images one by one.
    for f =  1:lic     % open f-th image-stack in the lif file.
%         if f == 3
%             thresForNuc =    .8      ;
%         end
        lifdat = bfopenOne(fToRead,f);
        %% Read metadata and Import image information from metadata
        [iminfo, Chinfo, chainfo] = read_and_import(lifdat,f);
        %% Detect a germline outline using ATS channel
        [lifdat, mip, mipBW, cyl_axis, cmask, pix] = detect_germline(lifdat, chainfo, Chinfo, iminfo, outline_channel);
        %% Detect Nuclei
        [nucs, thrNuc ] = detect_nuclei(lifdat, Chinfo, chainfo, iminfo, pix, sensi, f, lic, thresForNuc, mipBW, nrange);
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
                [blobrna, ~, mom1(g)] = DetectRNAIntron(thrNuc, rPla, pixr, iminfo, f, lic, g, thresForIntron);
                
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
                [DETrna, ~, mom2(g)] = DetectRNAExon(rPla, pixr, iminfo, f, lic, g, thresForExon);
                
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
        
        [ExATSrna_ATS, prna, ExATSrna] = validate(prna, ExATSrna, ChOrder, iminfo, nuc);
        
        
        
        
        
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