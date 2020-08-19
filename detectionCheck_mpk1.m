             %% Check detection for each channel

lic = size(af,1);
chname;
for f = 1:lic 
%         continue
%     end
    
    g=       2         ; % i-th channel to display


%     figure, imshow(af{f,8+g}*       5       )
    figure, imshow(af{f,9+g}*    15     ) %15 is constrast 
    hold on
    if ~isempty(af{f,g*2})
    plot(af{f,g*2}(:,1), af{f,g*2}(:,2), 'ro', 'markersize', 15, 'linewidth', 1);
    end
    
    
    fprintf('\n%d-th image', f);
      pause
    close all
end   
