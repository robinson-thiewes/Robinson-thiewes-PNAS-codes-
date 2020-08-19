function [] = VisualRNAdetect(af, g, marker, mag)
%%% This function visualizes detected RNA on top of original images. 
%%% The first input variable 'af' takes results data 'af' from smFISH
%%% analysis. The second input variable 'g' indicates the RNA channel to
%%% visualize (e.g. put 1 for g if you want to visualize 1st RNA channel).
%%% The third input is to determine how to show detected RNA dots. Use the
%%% same syntax as the function 'plot' (e.g. 'r+' is to show dots in red
%%% with crosses).

for f= 1:size(af,1)
    
%     if f == 4
%         continue
%     end

%         blankp = ones(size(af{f,8+g},1),size(af{f,8+g},2));
%         figure, imshow(blankp)
        figure, imshow(af{f,9+g}*mag)
        figure, imshow(af{f,9+g}*mag); hold on
        plot(af{f,g*2}(:,1), af{f,g*2}(:,2), marker, 'markersize', 10, 'linewidth', 2);
%         plot(prna{g}(:,1), prna{g}(:,2), 'ro', 'markersize', 20, 'linewidth', 3);
%         cirCiz = ceil(af{f,4}(:,6) / mean(af{f,4}(:,6)) * 15);
%         for i=1:length(af{2,4}(:,1))
%             plot(af{f,g*2}(i,1), af{f,g*2}(i,2), 'r.', 'markersize', cirCiz(i), 'linewidth', 3);
%         end
%         plot(af{f,g*2}(:,1), af{f,g*2}(:,2), 'm.', 'markersize', 23, 'linewidth', 3);

%%% ---- in case showing RNA IDs -------------
%     labes = 1:size(af{f,g*2},1);
%     labes = arrayfun(@num2str, labes, 'uni', 0);
%     text(af{f,g*2}(:,1), af{f,g*2}(:,2),labes, 'color', 'y', 'fontsize', 20)

%     plot(prna{g}(:,1), prna{g}(:,2), 'ro', 'markersize', 20, 'linewidth', 3);
%     labes = 1:length(prna{g}(:,1));
%     labes = arrayfun(@num2str, labes, 'uni', 0);
%     text(prna{g}(:,1), prna{g}(:,2),labes, 'color', 'y', 'fontsize', 20)
    fprintf('\n%d-th image',f); 
    drawnow;
    pause();
    close all
end