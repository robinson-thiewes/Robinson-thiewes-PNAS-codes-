function [lifdat, mip, mipBW, cyl_axis, cmask, pix] = detect_germline(lifdat, chainfo, Chinfo, iminfo, outline_channel)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Seperate planes
        % The last number on line below determines which color channel to use.
        PhaPla = lifdat{1,1}(chainfo(:) == find(Chinfo(:,2) == outline_channel), :);
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
        [cyl_axis, ~,cmask] = DefOutlineATS(mip, pix, iminfo);
        
        %%% align germline orientation (distal end to the left)
        if mean(cyl_axis(1:10,4)) > mean(cyl_axis(end-9:end,4))
            cyl_axis(:,1) = iminfo(2)-cyl_axis(:,1)+1;
            cmask = flip(cmask, 2);
            %lifdat is {1x4}, {1,1} is image pixels: {1,1}{:,1}, image name: {1,1}{:,2}
            lifdat{1,1}(:,1) = cellfun(@(x) flip(x,2), lifdat{1,1}(:,1), 'uni', 0); %flip images if germline is flipped
        end
end