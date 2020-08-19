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

%%% FPC{channel,2} data is a [#images (currImage) by 2] cell. {currImage,1} is a [#RNAs (currRNA) by 2]
%%% matrix , where {currRNA,1} is the RNA number or ID (if exists) and {currRNA,2} is
%%% the manual call of that RNA (1=yes,positive;  2=no,false-positive; 0 is
%%% yet unchecked manually or skipped.) Images are denoted by the currImage
%%% variable.

%%% Example: FPC{g,2}{currImage,1}[currRNA,2] = 2 means that the currRNA-th RNA in the currImage-th
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


%%% change g manually to change the channel; change line 117 for splitprobe

%% assert that an af variable is in the workspace
%lic - number of separate images
try
    lic = size(af,1);
catch
    warning('There is no af variable to be read.');
    return;
end
g=       2 ; % i-th channel to display - change this to change the channel




%% initialize FPC if there is none
%FPC - datastructure to store results
if exist('FPC','var') == 0
    % load is a boolean value containing true if this is a load of previous
    % progress
    load = false;
    FPC = cell(length(ChOrder),3);
end
if isempty(FPC{g, 1})
    FPC{g,1} = ChOrder(g);
    FPC{g,2} = cell(lic,3);

    
    for currImage = 1:lic
        if ~isempty(af{currImage,2*g})
            numRNA = numel(af{currImage,2*g}(:,8));
            FPC{g,2}{currImage,1} = zeros(numRNA, 2);
            if ChOrder(g) == 2
                FPC{g,2}{currImage,1}(:,1) = af{currImage,g*2}(:,8);
            else
                FPC{g,2}{currImage,1}(:,1) = 1:numRNA;
            end
        end
    end
    currImage = 1;
    
else
    load = true;
    % assert that if FBC does exist, the user is in the process of checking
    % false positives, and so the current image 'currImage' and the current RNA 'currRNA'
    % do exist
    %    assert(exist('currRNA','var') == 1, 'currRNA does not exist');
    %    assert(exist('currImage','var')==1, 'currImage does not exist');
    
    % If a FPC file has been opened and not newly created, the currImage image number
    % may not exist. Create currImage does not exist and set to the first image with
    % unchecked RNA.
    
    % If all RNA have been checked, the image number and RNA
    % number to start must be changed.
    
    if (exist('currImage','var'))== 0
        for currImage = lic:-1:1
            if ~isempty(af{currImage,g*2})
                currRNA = find(FPC{g,2}{currImage,1}(:,2) == 0, 1);
                if(~isempty(currRNA))
                    contImage = currImage;
                end
            end
        end
        if exist('contImage','var')== 1
            currImage = contImage;
        else
            currImage = 1;
        end
    end
end


disp('1 or y for yes, 2 or for no, 0 or s for skip, m for menu');
while currImage <= lic
    if ~isempty(af{currImage,2})
        %% show image
        numRNA = numel(FPC{g,2}{currImage,1}(:,1));
        fprintf('\n%d-th image', currImage);
        %     figure, imshow(af{currImage,8+g}*       5       )
        % change to 12+g for the splitprobe dataset
        figure, imshow(af{currImage,9+g}*    10      )
        hold on;
        for currRNAF = 1:numRNA
            if FPC{g,2}{currImage,1}(currRNAF,2)~=0
                plot(af{currImage,g*2}(currRNAF,1), af{currImage,g*2}(currRNAF,2), 'ro', 'markersize', 2, 'linewidth', 1);
            end
        end
        %%%%%%%figure, imshow(af{currImage,12+g}*    10      )
        %     hold on if ~isempty(af{currImage,g*2}) plot(af{currImage,g*2}(:,1),
        %     af{currImage,g*2}(:,2), 'ro', 'markersize', 15, 'linewidth', 1); end
        
        % find next RNA location
        fprintf('\n%d RNA locations to be shown \n', numRNA);
        if ~isempty(af{currImage,g*2})
            currRNA = find(FPC{g,2}{currImage,1}(:,2) == 0, 1);
            if(isempty(currRNA))
                currRNA = 1;
            end
            while currRNA <= numRNA
                hold on;
                %% add circles and ID numbers to the plot
                data = plot(af{currImage,g*2}(currRNA,1), af{currImage,g*2}(currRNA,2), 'ro', 'markersize', 15, 'linewidth', 1);
                %af(:,2) is prna prna(:,8) is RNA id
                num1 = text(af{currImage,g*2}(currRNA,1), af{currImage,g*2}(currRNA,2), num2str(FPC{g,2}{currImage,1}(currRNA,1)), 'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment','right');
                
                %% prompt for user input and act accordingly
                rnaidask = sprintf('Is RNA ID %d a true call? ', FPC{g,2}{currImage,1}(currRNA,1));
                in = input(rnaidask, 's');
                keepdata = false;
                switch in
                    case {'1','y', 'yes'}
                        FPC{g,2}{currImage,1}(currRNA,2) = 1;
                        keepdata = true;
                    case {'2','n', 'no'}
                        FPC{g,2}{currImage,1}(currRNA,2) = 2;
                        %                  case {'0','s', 'skip'}
                        keepdata = true;
                    case {'0','u', 'unsure'}
                        %s or save saves to a user prompted path
                    case{'s','save'}
                        %calculate ratios

                        for currImageF = 1:lic
                            if ~isempty(FPC{g,2}{currImageF,1})
                                % calculate columns 2 of FPC{g,2}
                                %# of correct calls
                                correct = sum(FPC{g,2}{currImageF,1}(:,2)==1);
                                %# of false positives
                                fp = sum(FPC{g,2}{currImageF,1}(:,2)==2);
                                %percentage of correct calls / determinable calls
                                FPC{g,2}{currImageF,2} = correct/(fp+correct);
                            end
                        end
                        %change path and save file
                        [FPCfile,FPCpath] = uiputfile('FPC.mat');
                        cd(FPCpath);
                        save(FPCfile,'FPC');
                        fprintf('\n---Saved FPC named %s to %s.---\n',FPCfile,FPCpath);
                        currRNA = currRNA-1;
                    case {'q','quit'}
                        % calculate ratios
                  
                        for currImageF = 1:lic
                            if ~isempty(FPC{g,2}{currImageF,1})
                                % calculate columns 2 of FPC{g,2}
                                %# of correct calls
                                correct = sum(FPC{g,2}{currImageF,1}(:,2)==1);
                                %# of false positives
                                fp = sum(FPC{g,2}{currImageF,1}(:,2)==2);
                                %percentage of correct calls / determinable calls
                                FPC{g,2}{currImageF,2} = correct/(fp+correct);
%                                 if sum(FPC{g,2}{currImageF,1}(:,2)==0) == 0
%                                     close all;
%                                     figure, imshow(af{currImage,9+g}*    10      )
%                                     hold on;
%                                     for currRNAF = 1:numRNA
%                                         if FPC{g,2}{currImageF,1}(currRNAF,2)~=0
%                                             plot(af{currImageF,g*2}(currRNAF,1), af{currImageF,g*2}(currRNAF,2), 'ro', 'markersize', 2, 'linewidth', 1);
%                                         end
%                                     end
%                                     uncalledAsk = sprintf('How many spots are not called? ');
%                                     uncalledSpotNum = input(uncalledAsk, 's');
%                                     FPC{g,2}{currImageF,3} = uncalledSpotNum;
%                                     
%                                 end
                            end
                        end
                        % quit script
                        return;
                    case {'back', 'b'}
                        delete(data);
                        currRNA = currRNA-2;
                    case {'menu', 'm'}
                        
                        
                        fprintf('\n%d-th image of %d.     %d-th of %d RNA locations to be shown.\n', ...
                            currImage, lic, currRNA, numRNA);
                        disp('1 or y for yes, 2 or n for no, 3 or u for unsure');
                        disp('back or b for backwards 1 RNA');
                        %                       disp('skip or s for skip');
                        disp('menu or m for menu');
                        disp('changeID or cID to change to a new RNA ID');
                        disp('changeI or cI to change to a new image');
                        disp('to change channels, please quit and manually change g');
                        currRNA = currRNA-1;
                        
                    case {'changeID', 'cID'}
                        fprintf('\n%d-th of %d RNA locations to be shown.\n', ...
                            currRNA, numRNA);
                        changeToVal = input('Change to which ID? ');
                        changeToRow = find(FPC{g,2}{currImage,1}(:,1) == changeToVal);
                        if isempty(changeToRow)
                            fprintf('That RNA ID could not be found.\n');
                        else
                            currRNA = changeToRow;
                        end
                        currRNA = currRNA-1;
                       
                    case {'changeImage', 'cI'}
                        fprintf('\n%d-th image of %d.\n', currImage, lic);
                        changeToImage = input('Change Image to which image? ');
                        if changeToImage > lic || changeToImage < 1
                            disp('That image could not be found.');
                            currRNA = currRNA-1;
                        else
                            currImage = changeToImage - 1;
                            currRNA = numRNA + 1;
                        end
                    otherwise
                        currRNA = currRNA-1;
                        disp('Unexpected answer. Please try again or press m to see the menu')
                end
                
                
                %% preparing for next RNA
                delete(data);
                if(keepdata)
                data = plot(af{currImage,g*2}(currRNA,1), af{currImage,g*2}(currRNA,2), 'ro', 'markersize', 2, 'linewidth', 1);
                end
                delete(num1);
                currRNA = currRNA + 1;
            end
        end
    end
    
    currImage = currImage+1;
    
    close all
end
%calculate ratios
for currImageF = 1:lic
    if ~isempty(FPC{g,2}{currImageF,1})
        % calculate columns 2 of FPC{g,2}
        %# of correct calls
        correct = sum(FPC{g,2}{currImageF,1}(:,2)==1);
        %# of false positives
        fp = sum(FPC{g,2}{currImageF,1}(:,2)==2);
        %percentage of correct calls / determinable calls
        FPC{g,2}{currImageF,2} = correct/(fp+correct);
    end
end

function FPC = calculate_ratios(FPC, lic, g, currImageF)
for currImageF = 1:lic
    if ~isempty(FPC{g,2}{currImageF,1})
        % calculate columns 2 of FPC{g,2}
        %# of correct calls
        correct = sum(FPC{g,2}{currImageF,1}(:,2)==1);
        %# of false positives
        fp = sum(FPC{g,2}{currImageF,1}(:,2)==2);
        %percentage of correct calls / determinable calls
        FPC{g,2}{currImageF,2} = correct/(fp+correct);
    end
end
end

% % function FPC = calcRatios(FPC,lic)
% % for currImageF = 1:lic
% %     if ~isempty(FPC{g,2}{currImageF,1})
% %         % calculate columns 2 of FPC{g,2}
% %         %# of correct calls
% %         correct = sum(FPC{g,2}{currImageF,1}(:,2)==1);
% %         %# of false positives
% %         fp = sum(FPC{g,2}{currImageF,1}(:,2)==2);
% %         %percentage of correct calls / determinable calls
% %         FPC{g,2}{currImageF,2} = correct/(fp+correct);
% %     end
% % end
% % end


%{ %% old code
% %% return to skipped values
% fprintf('\n\n     Displaying skipped values \n\n');
% for currImage = 1:lic
%    while(~all(FPC{g,2}{currImage,1}(:,2)))
%        fprintf('Displaying image %d.\n', currImage);
%        currRNA = find(FPC{g,2}{currImage,1}(:,2) == 0, 1);
%        DispAsk(af, currImage, g, currRNA, FPC{g,2});
%
%    end
% end

% %% function to display RNAs and ask for user input
% function [currImage, currRNA, FPC{g,2}] = DispAsk(af, currImage, g, currRNA, FPC{g,2})
%             %% add circles and ID numbers
%             hold on
%             data = plot(af{currImage,g*2}(currRNA,1), af{currImage,g*2}(currRNA,2), 'ro', 'markersize', 15, 'linewidth', 1);
%             %af(:,2) is prna prna(:,8) is RNA id
%             num1 = text(af{currImage,g*2}(currRNA,1), af{currImage,g*2}(currRNA,2), num2str(FPC{g,2}{currImage,1}(currRNA,1)), 'VerticalAlignment', 'bottom', ...
%                 'HorizontalAlignment','right');
%
%             %% prompt for user input and act accordingly
%             rnaidask = sprintf('\nIs RNA ID %d a true call?\n', FPC{g,2}{currImage,1}(currRNA,1));
%             in = input(rnaidask, 's');
%             switch in
%                 case {'1','y', 'yes'}
%                     FPC{g,2}{currImage,1}(currRNA,2) = 1;
%                 case {'2','n', 'no'}
%                     FPC{g,2}{currImage,1}(currRNA,2) = 2;
%                 case {'0','s', 'skip'}
%                 case {'back', 'b'}
%                     currRNA = currRNA-2;
%                 case {'menu', 'm'}
%
%                     fprintf('\n%d-th image of %d.     %d-th of %d RNA locations to be shown.\n', ...
%                         currImage, lic, currRNA, numRNA);
%                     disp('1 or y for yes, 2 or for no, 3 or u for unsure');
%                     disp('back or b for backwards 1 RNA');
%                     disp('skip or s for skip');
%                     disp('menu or m for menu');
%                     disp('changeID or cID to change to a new RNA ID');
% %                   disp('changeImage or cI to change to a new image');
% %                     %this is not functional
%                     currRNA = currRNA-1;
%                 case {'changeID', 'cID'}
%                     fprintf('\n%d-th of %d RNA locations to be shown.\n', ...
%                         currRNA, numRNA);
%                     changeToVal = input('Change to which ID?');
%                     changeToRow = find(FPC{g,2}{currImage,1}(:,1) == changeToVal);
%                     if isempty(changeToRow)
%                         disp('That RNA ID could not be found.');
%                     else
%                         currRNA = changeToRow;
%                     end
%                     currRNA = currRNA-1;
%
% %                 %this is not yet functional
% %                 case {'changeImage', 'cI'}
% %                     fprintf('\n%d-th image of %d.\n', currImage, lic);
% %                     changeToImage = input('Change Image to which
% %                     image?'); if changeToImage > lic
% %                         disp('That image could not be found.');
% %                     else
% %                         currImage = changeToImage; currRNA = numRNA + 1;
% %                     end
%                     currRNA = currRNA-1;
%                 otherwise
%                     currRNA = currRNA-1;
%                     warning('Unexpected answer. Please try again or press m to see the menu')
%             end
%
%             %% preparing for next RNA
%
%             hold off
%             delete(data);
%             delete(num1);
%             currRNA = currRNA + 1;
% end

% FPC = FPC{g,2}
%}

