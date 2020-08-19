function [perATSpCR, meanATSpCR, stdATSpCR]=PerCellATS(aff, ATSch, iminfo)
% [perATSpCR, meanATSpCR, stdATSpCR, perATSpCR_exon ,meanATSpCR_exon, stdATSpCR_exon]=PerCellATS(aff, ATSch, iminfo)
%%% This code draws plot for % of cells with ATS against distance (um) from distal end.

rc = ATSch; %ATS channel
% rc = 2;
rc = rc*2-1;
ic = ATSch*2;
crs=5;
xcr = 0:crs:60;
bins = length(xcr);
xc = ATSch*4;
xn = xc-1;
perATSpCR = zeros(size(aff,1), bins); %take nuclear information
perATSpCR_exon = zeros(size(aff,1), bins);

for i=1:size(aff,1) %i = # of images
%     [~,ia,ib] = intersect(aff{i,ic}(:,4), aff{i,xn}(:,8), 'rows', 'stable'); %find intersection between nuc ID listed in intron ATS info and exon nuc info
%     exon_ATS = aff{i,ic}(ia,:); %index interesected nucIDs
%     exon_ATS = aff{i,ic}(:,:);
%     aff{i,rc}(:,10) = aff{i,xn}(:,8);
%     exon_ATS_only = exon_ATS(:,10)==0;
%     exon_ATS(exon_ATS_only,:) = []; %makes list of exons matching intron detected in intron ATS list
%     counts = countmember(aff{i,rc}(:,10), exon_ATS(:,4)); %counts number of exon ATS per nucID
%     nucs = aff{i,rc}(:,:);
%     nucs(:,11) = counts;
    for j=1:bins     % 1 - nth cell row
        xr = [crs*(j-1)/iminfo(6) crs*j/iminfo(6)];
        allnucs = sum(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2));
        if allnucs == 0
            allnucs = 1;
        end
        perATSpCR(i,j) = sum(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2) & aff{i,rc}(:,9) > 0)  / ...
            allnucs * 100;
%         perATSpCR_exon(i,j)=sum(nucs(:,1) >= xr(1) & nucs(:,1) < xr(2) & nucs(:,11) > 0)  / ...
%             allnucs * 100;
    end
end

meanATSpCR = mean(perATSpCR);
stdATSpCR = std(perATSpCR);
% meanATSpCR_exon = mean(perATSpCR_exon);
% stdATSpCR_exon = std(perATSpCR_exon);
% ------ bar graph with multiple strains
% figure
% hold on
% 
% bar(xcr,[meanATSpCR_N2', meanATSpCR_q46', meanATSpCR_q385', meanATSpCR_q411']);
% errorbar([bk'-1.35, bk'-.45, bk'+.45, bk'+1.35], [meanATSpCR_N2', meanATSpCR_q46', meanATSpCR_q385', meanATSpCR_q411'], ...
%     [stdATSpCR_N2'/sqrt(78), stdATSpCR_q46'/sqrt(20), stdATSpCR_q385'/sqrt(20), stdATSpCR_q411'/sqrt(20)], 'k.', 'linewidth', 2);
% axis([ -4 64 0 80 ])



%----------bargraph
figure('pos', [300 200 350 500])
hold on
bar(xcr+crs/2, meanATSpCR)
errorbar(xcr+crs/2, meanATSpCR, stdATSpCR/sqrt(size(perATSpCR,1)), 'k.', 'linewidth', 2);
xticks(0:5:60)
axis([ 0 60 0 100])
axisH = axis;
%[nAxis] = addTopXAxis(axisH,properties...)
axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
% box on
ylabel('% cell with \itmpk-1\rm ATS by intron probe' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);
% hold off
% 
% 
% %% graph of exon detected ATS
% figure('pos', [300 200 350 500])
% 
% hold off
% h1=bar(xcr+crs/2, meanATSpCR_exon);
% h1.FaceColor = [1 0 .5];
% hold on
% errorbar(xcr+crs/2, meanATSpCR_exon, stdATSpCR_exon/sqrt(size(perATSpCR_exon,1)), 'k.', 'linewidth', 2);
% axis([ -2 60 0 100])
% xticks(0:5:60)
% axisH = axis;
% %[nAxis] = addTopXAxis(axisH,properties...)
% axisH = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
% box on
% ylabel('% \itmpk-1\rm ATS by exon' , 'fontsize',15);
% xlabel('\itum\rm from distal end', 'fontsize',15);


% figure('pos', [300 200 350 500])
% hold on
% bar(xcr+crs/2, meanATSpCR_exon)
% errorbar(xcr+crs/2, meanATSpCR_exon, stdATSpCR_exon/sqrt(size(perATSpCR_exon,1)), 'k.', 'linewidth', 2);
% xticks(0:5:60)
% axis([ 0 60 0 100])
% axisI = axis;
% %[nAxis] = addTopXAxis(axisH,properties...)
% axisI = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
% % box on
% ylabel('% cell with \itmpk-1\rm ATS by exon probe' , 'fontsize',15);
% xlabel('\itum\rm from distal end', 'fontsize',15);
% hold off

%% graph of intron and exon detected ATS

% figure('pos', [300 200 350 500])
% hold off
% hbar = bar(xcr+crs/2, [meanATSpCR',meanATSpCR_exon']);
% hold on
% pause(0.1)
% errorbar(xcr+crs/2+hbar(1).XOffset, meanATSpCR',stdATSpCR'/sqrt(size(perATSpCR,1)), 'k.', 'linewidth',2);  
% errorbar(xcr+crs/2+hbar(2).XOffset, meanATSpCR_exon',stdATSpCR_exon'/sqrt(size(perATSpCR_exon,1)), 'k.', 'linewidth',2); 
% xticks(0:5:60);
% axis([ -2 60 0 100 ]);
% axisJ = axis;
% %[nAxis] = addTopXAxis(axisH,properties...)
% axisJ = addTopXAxis('expression', '(argu/4.4)', 'xLabStr', 'gcd in window');
% % legend('Intron detected ATS', 'Exon detected ATS', 'location', 'northeast')
% ylabel('% cell with \itmpk-1b\rm ATS' , 'fontsize',15);
% xlabel('\itum\rm from distal end', 'fontsize',15);

