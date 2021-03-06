function PerCellATScirc(aff, ATSch, iminfo, limx)
%%% This code draws plot for % of cells with ATS against distance (um) from distal end.

rc = ATSch;
rc = rc*2-1;

crs=3;
xcr = 0:crs:25;
bins = length(xcr);

perATSpCR1 = zeros(size(aff,1), bins-1);
perATSpCR2 = zeros(size(aff,1), bins-1);

for i=1:size(aff,1)
    mn = ( max(aff{i,1}(:,3))+min(aff{i,1}(:,3)) )/2;
    
    
    for j=1:bins-1     % 1 - 7th cell row        
        xr = [xcr(j)/iminfo(6) xcr(j+1)/iminfo(6)];
        
        allnucs = sum(aff{i,rc}(:,2) >= xr(1) & aff{i,rc}(:,2) < xr(2) & aff{i,rc}(:,1) < limx/iminfo(6) & aff{i,rc}(:,3) < mn);
        if allnucs == 0
            allnucs = 1;
        end
        
        perATSpCR1(i,j) = sum(aff{i,rc}(:,2) >= xr(1) & aff{i,rc}(:,2) < xr(2) & aff{i,rc}(:,1) ...
            < limx/iminfo(6) & aff{i,rc}(:,3) < mn & aff{i,rc}(:,8) > 0) / allnucs * 100;
    end
    
    for j=1:bins-1     % 1 - 7th cell row
        xr2 = [xcr(bins-j)/iminfo(6) xcr(bins-j+1)/iminfo(6)];
        allnucs = sum(aff{i,rc}(:,2) >= xr2(1) & aff{i,rc}(:,2) < xr2(2) & aff{i,rc}(:,1) < limx/iminfo(6) & aff{i,rc}(:,3) >= mn);
        
        if allnucs == 0
            allnucs = 1;
        end
        
        perATSpCR2(i,j) = sum(aff{i,rc}(:,2) >= xr2(1) & aff{i,rc}(:,2) < xr2(2) & aff{i,rc}(:,1) ...
            < limx/iminfo(6) & aff{i,rc}(:,3) >= mn & aff{i,rc}(:,8) > 0) / allnucs * 100;
    end
end

perATSpCR = [perATSpCR1 perATSpCR2];


meanATSpCR = mean(perATSpCR);
stdATSpCR = std(perATSpCR);

%------ bar graph with multiple strains
% figure
% hold on
% 
% bar(xcr,[meanATSpCR_N2', meanATSpCR_q46', meanATSpCR_q385', meanATSpCR_q411']);
% errorbar([bk'-1.35, bk'-.45, bk'+.45, bk'+1.35], [meanATSpCR_N2', meanATSpCR_q46', meanATSpCR_q385', meanATSpCR_q411'], ...
%     [stdATSpCR_N2'/sqrt(78), stdATSpCR_q46'/sqrt(20), stdATSpCR_q385'/sqrt(20), stdATSpCR_q411'/sqrt(20)], 'k.', 'linewidth', 2);
% axis([ -4 64 0 80 ])

xcr(end)=[];
xcr = [xcr xcr+max(xcr)+crs];


%----------bargraph
figure('pos', [300 200 350 500])

hold off
bar(xcr+crs/2, meanATSpCR)
hold on
errorbar(xcr+crs/2, meanATSpCR, stdATSpCR/sqrt(size(perATSpCR,1)), 'k.', 'linewidth', 2);
xticks(0:10:100)
axis([ -2 62 0 100 ])
box on
ylabel('% cell with \itsygl-1\rm ATS (0-30um from distal end)' , 'fontsize',15);
xlabel('\itum\rm around the gonad', 'fontsize',15);


% 
% %%% viewing % cells w/ ATS (circumference) plot for individual germlines
% for i=1:size(perATSpCR,1)
%     figure
%     bar(1:size(perATSpCR,2),perATSpCR(i,:));
%     box on
%     pause
%     close all
% end
% 
% 










