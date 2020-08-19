 %% Detect Nuclei
%written by Changhwan, editted by JM
function [nucs, thrNuc ] = detect_nuclei(lifdat, Chinfo, chainfo, iminfo, pix, sensi, f, lic, thresForNuc, mipBW, nrange)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % seperate Nuc channel to detect nuclei.
        NucPla = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == 1), :);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Detect nuclei and remove one nucleus if two are overlapped.
        % NDNuc saves not overlapped Nuclei that are in DetNuc (detected Nuc).
        tnran = round(nrange ./ iminfo(5:6,1)');
        tol = round( 1.2 / mean(iminfo(5:6,1))); %%%%% tolerance for overlapped circles (um).
        
        % Detect nuclei from image slices and make nuclear sphere.
        
        
        [Nuc, thrNuc, aftNuc, ~] = DetectNucleus(NucPla, iminfo, pix, tnran, sensi, tol, f, lic, thresForNuc);
        
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
        
        