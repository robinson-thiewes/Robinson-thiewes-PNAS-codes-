%%% This code will display RNAs from a given af file individually and
%%% prompt the user to manually call whether they are indeed RNAs or just
%%% false positives.

%%% This was adapted from smFISH_analysis. The af and ChOrder variables from the final
%%% workspace of smFISH_v1_02.m etc. is required to run this code.

%%% Results are stored in FPC. (False positive check)
%%% FPC is a [#channels (g) by 3] cell. {g,1} is 
%%% the name of the channel. {g,2} contains channel RNA false
%%% postive check data. {g,3} is the percentage of RNA that was manually
%%% checked to call correctly over RNA overall manually checked. Channels
%%% are noted by the g variable.

%%% FPC{channel,2} data is a [#images (fj) by 2] cell. {fj,1} is a [#RNAs (j) by 2]
%%% matrix , where {j,1} is the RNA number or ID (if exists) and {j,2} is
%%% the manual call of that RNA (1=yes,positive;  2=no,false-positive; 0 is
%%% yet unchecked manually or skipped.) Images are denoted by the fj
%%% variable.

%%% Example: FPC{g,2}{fj,1}[j,2] = 2 means that the j-th RNA in the fj-th
%%% image in the g-th channel is not an actual RNA (has been falsely
%%% called).


%%% A variety of commands are available during use.
%%% 1 or y or yes for yes, 2 or n or no for no, 3 or u for unsure
%%% back or b for backwards 1 RNA
%%% s or save for save (will save the current FPC file to a prompted
%%% location)
%%% q or quit for quit
%%% menu or m for menu (shows available commands)
%%% changeID or cID to change to a new RNA ID
%%% changeImage or cI to change to a new image


%%% This code has not been thoroughly debugged yet; please contact John at
%%% jmccloskey2@wisc.edu with any questions or bugs.

%%% The long comment block in the end is for reworking possible areas of
%%% the code and can be deleted.

%% assert that an af variable is in the workspace
%lic - number of separate images
try
    lic = size(af,1);
catch
    warning('There is no af variable to be read.');
    return;
end
g=       1 ; % i-th channel to display - change this to change the channel


%% initialize FPC if there is none
%FPC - datastructure to store results
if exist('FPC','var') == 0
    FPC = cell(length(ChOrder),3);
    FPC{g,1} = ChOrder(g);
    FPC{g,2} = cell(lic,2);
    
    for jf = 1:lic
        if ~isempty(af{jf,2})
            numRNA = numel(af{jf,2}(:,8));
            FPC{g,2}{jf,1} = zeros(numRNA, 2);
            if ChOrder(g) == 2
                FPC{g,2}{jf,1}(:,1) = af{jf,g*2}(:,8);
            else
                FPC{g,2}{jf,1}(:,1) = 1:numRNA;
            end
        end
    end
    jf = 1;
    
else
    % assert that if FBC does exist, the user is in the process of checking
    % false positives, and so the current image 'jf' and the current RNA 'j'
    % do exist
%    assert(exist('j','var') == 1, 'j does not exist');
%    assert(exist('jf','var')==1, 'jf does not exist');

% If a FPC file has been opened and not newly created, the jf image number 
% may not exist. Create jf does not exist and set to the first image with
% unchecked RNA. 

% If all RNA has been checked, the image number and RNA
% number to start must be changed.
          
    if (exist('jf','var'))== 0
        for jf = lic:-1:1
            if ~isempty(af{jf,g*2})
                j = find(FPC{g,2}{jf,1}(:,2) == 0, 1);
                if(~isempty(j))
                    contImage = jf;
                end
            end
        end
        if exist('contImage','var')== 1
            jf = contImage;
        else
            jf = 1;
        end
    end
end


disp('1 or y for yes, 2 or for no, 0 or s for skip, m for menu');
while jf <= lic
    if ~isempty(af{jf,2})
        %% show image
        numRNA = numel(FPC{g,2}{jf,1}(:,1));
        fprintf('\n%d-th image', jf);
        %     figure, imshow(af{jf,8+g}*       5       )
        % change to 12+g for the splitprobe dataset
        figure, imshow(af{jf,9+g}*    10      )
        %%%%%%%figure, imshow(af{jf,12+g}*    10      )
        %     hold on if ~isempty(af{jf,g*2}) plot(af{jf,g*2}(:,1),
        %     af{jf,g*2}(:,2), 'ro', 'markersize', 15, 'linewidth', 1); end
        
        % find next RNA location
        fprintf('\n%d RNA locations to be shown \n', numRNA);
        if ~isempty(af{jf,g*2})
            j = find(FPC{g,2}{jf,1}(:,2) == 0, 1);
            if(isempty(j))
                j = 1;
            end
            while j <= numRNA
                %% add circles and ID numbers to the plot
                hold on
                data = plot(af{jf,g*2}(j,1), af{jf,g*2}(j,2), 'ro', 'markersize', 15, 'linewidth', 1);
                %af(:,2) is prna prna(:,8) is RNA id
                num1 = text(af{jf,g*2}(j,1), af{jf,g*2}(j,2), num2str(FPC{g,2}{jf,1}(j,1)), 'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment','right');
                
                %% prompt for user input and act accordingly
                rnaidask = sprintf('Is RNA ID %d a true call? ', FPC{g,2}{jf,1}(j,1));
                in = input(rnaidask, 's');
                switch in
                    case {'1','y', 'yes'}
                        FPC{g,2}{jf,1}(j,2) = 1;
                    case {'2','n', 'no'}
                        FPC{g,2}{jf,1}(j,2) = 2;
                        %                  case {'0','s', 'skip'}
                    case {'0','u', 'unsure'}
                    %s or save saves to a user prompted path
                    case{'s','save'}
                        %calculate ratios
                        for ff = 1:lic
                            if ~isempty(FPC{g,2}{ff,1})
                                % calculate columns 2 of FPC{g,2}
                                %# of correct calls
                                correct = sum(FPC{g,2}{ff,1}(:,2)==1);
                                %# of false positives
                                fp = sum(FPC{g,2}{ff,1}(:,2)==2);
                                %percentage of correct calls / determinable calls
                                FPC{g,2}{ff,2} = correct/(fp+correct);
                            end
                        end
                        %change path and save file
                        [FPCfile,FPCpath] = uiputfile('FPC.mat');
                        cd(FPCpath);
                        save(FPCfile,'FPC');
                        fprintf('\n---Saved FPC named %s to %s.---\n',FPCfile,FPCpath);
                        j = j-1;
                    case {'q','quit'}
                        % calculate ratios
                        for ff = 1:lic
                            if ~isempty(FPC{g,2}{ff,1})
                                % calculate columns 2 of FPC{g,2}
                                %# of correct calls
                                correct = sum(FPC{g,2}{ff,1}(:,2)==1);
                                %# of false positives
                                fp = sum(FPC{g,2}{ff,1}(:,2)==2);
                                %percentage of correct calls / determinable calls
                                FPC{g,2}{ff,2} = correct/(fp+correct);
                            end
                        end
                        % quit script
                        return;
                    case {'back', 'b'}
                        j = j-2;
                    case {'menu', 'm'}
                        
                        fprintf('\n%d-th image of %d.     %d-th of %d RNA locations to be shown.\n', ...
                            jf, lic, j, numRNA);
                        disp('1 or y for yes, 2 or n for no, 3 or u for unsure');
                        disp('back or b for backwards 1 RNA');
 %                       disp('skip or s for skip');
                        disp('menu or m for menu');
                        disp('changeID or cID to change to a new RNA ID');
                        disp('changeI or cI to change to a new image');
                        disp('to change channels, please quit and manually change g');
                        j = j-1;
                    case {'changeID', 'cID'}
                        fprintf('\n%d-th of %d RNA locations to be shown.\n', ...
                            j, numRNA);
                        changeToVal = input('Change to which ID? ');
                        changeToRow = find(FPC{g,2}{jf,1}(:,1) == changeToVal);
                        if isempty(changeToRow)
                            fprintf('That RNA ID could not be found.\n');
                        else
                            j = changeToRow;
                        end
                        j = j-1;
                    case {'changeImage', 'cI'}
                        fprintf('\n%d-th image of %d.\n', jf, lic);
                        changeToImage = input('Change Image to which image? ');
                        if changeToImage > lic || changeToImage < 1
                            disp('That image could not be found.');
                            j = j-1;
                        else
                            jf = changeToImage - 1;
                            j = numRNA + 1;
                        end
                    otherwise
                        j = j-1;
                        disp('Unexpected answer. Please try again or press m to see the menu')
                end
                
                
                %% preparing for next RNA
                
                hold off
                delete(data);
                delete(num1);
                j = j + 1;
            end
        end
    end
    
    jf = jf+1;
    
    close all
end
                        %calculate ratios
                        for ff = 1:lic
                            if ~isempty(FPC{g,2}{ff,1})
                                % calculate columns 2 of FPC{g,2}
                                %# of correct calls
                                correct = sum(FPC{g,2}{ff,1}(:,2)==1);
                                %# of false positives
                                fp = sum(FPC{g,2}{ff,1}(:,2)==2);
                                %percentage of correct calls / determinable calls
                                FPC{g,2}{ff,2} = correct/(fp+correct);
                            end
                        end


% % function FPC = calcRatios(FPC,lic)
% % for ff = 1:lic
% %     if ~isempty(FPC{g,2}{ff,1})
% %         % calculate columns 2 of FPC{g,2}
% %         %# of correct calls
% %         correct = sum(FPC{g,2}{ff,1}(:,2)==1);
% %         %# of false positives
% %         fp = sum(FPC{g,2}{ff,1}(:,2)==2);
% %         %percentage of correct calls / determinable calls
% %         FPC{g,2}{ff,2} = correct/(fp+correct);
% %     end
% % end
% % end


%{ %% old code
% %% return to skipped values
% fprintf('\n\n     Displaying skipped values \n\n');
% for jf = 1:lic
%    while(~all(FPC{g,2}{jf,1}(:,2)))
%        fprintf('Displaying image %d.\n', jf);
%        j = find(FPC{g,2}{jf,1}(:,2) == 0, 1);
%        DispAsk(af, jf, g, j, FPC{g,2});
%
%    end
% end

% %% function to display RNAs and ask for user input
% function [jf, j, FPC{g,2}] = DispAsk(af, jf, g, j, FPC{g,2})
%             %% add circles and ID numbers
%             hold on
%             data = plot(af{jf,g*2}(j,1), af{jf,g*2}(j,2), 'ro', 'markersize', 15, 'linewidth', 1);
%             %af(:,2) is prna prna(:,8) is RNA id
%             num1 = text(af{jf,g*2}(j,1), af{jf,g*2}(j,2), num2str(FPC{g,2}{jf,1}(j,1)), 'VerticalAlignment', 'bottom', ...
%                 'HorizontalAlignment','right');
%
%             %% prompt for user input and act accordingly
%             rnaidask = sprintf('\nIs RNA ID %d a true call?\n', FPC{g,2}{jf,1}(j,1));
%             in = input(rnaidask, 's');
%             switch in
%                 case {'1','y', 'yes'}
%                     FPC{g,2}{jf,1}(j,2) = 1;
%                 case {'2','n', 'no'}
%                     FPC{g,2}{jf,1}(j,2) = 2;
%                 case {'0','s', 'skip'}
%                 case {'back', 'b'}
%                     j = j-2;
%                 case {'menu', 'm'}
%
%                     fprintf('\n%d-th image of %d.     %d-th of %d RNA locations to be shown.\n', ...
%                         jf, lic, j, numRNA);
%                     disp('1 or y for yes, 2 or for no, 3 or u for unsure');
%                     disp('back or b for backwards 1 RNA');
%                     disp('skip or s for skip');
%                     disp('menu or m for menu');
%                     disp('changeID or cID to change to a new RNA ID');
% %                   disp('changeImage or cI to change to a new image');
% %                     %this is not functional
%                     j = j-1;
%                 case {'changeID', 'cID'}
%                     fprintf('\n%d-th of %d RNA locations to be shown.\n', ...
%                         j, numRNA);
%                     changeToVal = input('Change to which ID?');
%                     changeToRow = find(FPC{g,2}{jf,1}(:,1) == changeToVal);
%                     if isempty(changeToRow)
%                         disp('That RNA ID could not be found.');
%                     else
%                         j = changeToRow;
%                     end
%                     j = j-1;
%
% %                 %this is not yet functional
% %                 case {'changeImage', 'cI'}
% %                     fprintf('\n%d-th image of %d.\n', jf, lic);
% %                     changeToImage = input('Change Image to which
% %                     image?'); if changeToImage > lic
% %                         disp('That image could not be found.');
% %                     else
% %                         jf = changeToImage; j = numRNA + 1;
% %                     end
%                     j = j-1;
%                 otherwise
%                     j = j-1;
%                     warning('Unexpected answer. Please try again or press m to see the menu')
%             end
%
%             %% preparing for next RNA
%
%             hold off
%             delete(data);
%             delete(num1);
%             j = j + 1;
% end

% FPC = FPC{g,2}
%}

