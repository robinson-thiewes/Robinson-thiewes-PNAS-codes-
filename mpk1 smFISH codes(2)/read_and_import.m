%% Read metadata and Import image information from metadata
%from smFISH_v1_02_mpk1_v1 by Sarah and Changhwan, editted by John
function [lifdat, iminfo, Chinfo, chainfo] = read_and_import(fToRead,f)

    lifdat = bfopenOne(fToRead,f);
    %%% read metadata
    meta = lifdat{1, 4};
    Nch = meta.getPixelsSizeC(f-1).getValue();
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
    %%% Import image information from metadata
    % find Nuc, phase or trx channel by color (channel label does not work).
    % Chinfo: |col 1: color info in the order of acquisition|
    
    %%Import image information from metadata
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
end
