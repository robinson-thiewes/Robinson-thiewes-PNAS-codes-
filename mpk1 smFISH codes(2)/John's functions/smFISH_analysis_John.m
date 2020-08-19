lic = size(af,1);

for f = 1:lic
    
%     if f == 8
%         continue
%     end
    
    g=       1 ; % i-th channel to display
    numRNA = numel(af{1,2}(:,8));
callTF = zeros(numRNA, 2);
callTF(:,1) = af{f,g*2}(:,8);

%     figure, imshow(af{f,8+g}*       5       )
    figure, imshow(af{f,9+g}*    10      )
%     hold on
%     if ~isempty(af{f,g*2})
%     plot(af{f,g*2}(:,1), af{f,g*2}(:,2), 'ro', 'markersize', 15, 'linewidth', 1);
%     end

    fprintf('\n%d RNA locations to be shown', numRNA);
    if ~isempty(af{f,g*2})
        j = 1;
        while j < numRNA
        
        hold on
        data = plot(af{f,g*2}(j,1), af{f,g*2}(j,2), 'ro', 'markersize', 15, 'linewidth', 1); 
        %af(:,2) is prna
        %prna(:,8) is RNA id
        num1 = text(af{f,g*2}(j,1), af{f,g*2}(j,2), num2str(af{f,g*2}(j,8)), 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment','right');
        rnaidask = sprintf('\nIs RNA ID %d a true call?\n', af{f,g*2}(j,8));
        in = input(rnaidask, 's');
        switch in
            case {'1','y', 'yes'}
                callTF(j,2) = 1;
            case {'2','n', 'no'}
                callTF(j,2) = 2;
            case {'3','u', 'unsure'}
                callTF(j,2) = 3;
            case {'back', 'b'}
                j = j-2;
            case {'skip', 's'}
                j = j;
            case {'menu', 'm'}
                j = j-1;
                 disp('1 or y for yes, 2 or for no, 3 or u for unsure');
                 disp('back or b for backwards 1 RNA');
                 disp('skip or s for skip');
                 disp('menu or m for menu');
                 disp('change or c to change to a new RNA ID');
            case {'change', 'c'}
                try
                changeToVal = input('Change to which value?');
                changeToRow = find(af{f,g*2}(:,8) == changeToVal);
                catch
                    warning('That RNA ID could not be found.');
                end
                j = changeToRow - 1;
            otherwise
                j = j-1;
                warning('Unexpected answer. Please try again or press m to see the menu')
        end

                
        
        hold off
        delete(data);
        delete(num1);
        j = j + 1;
        end
    end

    fprintf('\n%d-th image', f);
    pause
    close all
end   