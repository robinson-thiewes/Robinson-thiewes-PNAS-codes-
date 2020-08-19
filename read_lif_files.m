    %% Read lif files
    %by Sarah and CHW, editted by John to adapt to smFISH_v1_02_mpk1_v1_2
function [nbf, liflist] = read_lif_files(MASTER_Image_path)
    

    %%% lif file format: separate diff condition by ' ' (space). Order or # of
    % statements does not matter. e.g. 040914 rrf-1 emptyRNAi_na lag-3_na.lif.
    % liflist = row1: lif file name w/o '.lif'; row2: name with sorted; conditions; row3: strain numbering.
    [liflist] = dir(strcat(MASTER_Image_path,'*.lif'));
    liflist = liflist(~[liflist.isdir]);
    liflist = {liflist.name}';
    nbf = numel(liflist);

end