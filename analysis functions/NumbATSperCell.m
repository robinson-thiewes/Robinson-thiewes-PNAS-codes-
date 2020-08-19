function NumbATSperCell(aff, ATSch, y_max,iminfo)




rc = ATSch;
rc = rc*2-1;

crs=5;
xcr = 0:crs:60;
bins = length(xcr);

perATSpCR = zeros(size(aff,1), 13);
for i=1:size(aff,1)
    for j=1:bins     % 1 - 13th cell row
        
        chooser = 0;   %%% 0 when considering all nucs, 1 when considering only nucs with any # ATS.
        
        xr = [crs*(j-1)/iminfo(6) crs*j/iminfo(6)];
        allnucs = sum(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2) & aff{i,rc}(:,8) > chooser - 1);
        if allnucs == 0
            allnucs = 1;
        end
        perATSpCR(i,j) = sum(aff{i,rc}(aff{i,rc}(:,1) >= xr(1) & aff{i,rc}(:,1) < xr(2) ...
            & aff{i,rc}(:,8) > chooser - 1, 8)) / allnucs;
        
    end
end

meanATSpCR = mean(perATSpCR);
stdATSpCR = std(perATSpCR);


%----------bargraph
figure('pos', [300 200 350 500])

hold off
bar(xcr+crs/2, meanATSpCR)
hold on
errorbar(xcr+crs/2, meanATSpCR, stdATSpCR/sqrt(size(perATSpCR,1)), 'k.', 'linewidth', 2);
axis([ -4 64 0 y_max ])
xticks(0:10:100)
box on
ylabel('# \itsygl-1\rm ATS per cell' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);






