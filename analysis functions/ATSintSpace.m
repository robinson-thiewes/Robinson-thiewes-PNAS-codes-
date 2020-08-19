function ATSintSpace(aff, ATSch, iminfo, y_max, atsNucAll)




rc = ATSch;
rc = rc*2-1;
divider = 13;     %%% 13 for old sygl-1 N2, bin number
ATSdistDE = cell(size(aff,1),1);
for i=1:size(aff,1)
    ATSp = aff{i,rc}(aff{i,rc}(:,8) > 0 & aff{i,rc}(:,9) > 0,:);% & aff{i,rc}(:,10) > 0,:);
%         xbor = floor(lengthChainBorder(i,2)/iminfo(6) *      0.9    );
    if aff{i,5}(1,1) > aff{i,5}(end,1)
        idx = size(aff{i,5},1);
    else
        idx = 1;
    end
    
    DTC = [aff{i,5}(idx,1) aff{i,5}(idx,2) aff{i,5}(idx,3)*iminfo(7)/iminfo(6)]; 
    [~, dd] = knnsearch(DTC, ATSp(:,1:3), 'k', 1);
    
    ATSdistDE{i,1}(:,1) = dd;
    ATSdistDE{i,1}(:,2) = ATSp(:,9);
    ATSdistDE{i,1}(:,4) = ATSp(:,9)/divider;
end

ATSdistDEmat = cell2mat(ATSdistDE);


%gets mean of transcriptional activity
gap = 1;
mLine = zeros(1,30);
seL = zeros(1,30);
for i=1:30
    cpool = ATSdistDEmat(ATSdistDEmat(:,1) > (i-1)/iminfo(6)*2 & ATSdistDEmat(:,1) < i/iminfo(6)*2,2);
    cpool = atsNucAll(atsNucAll(:,1) > (i-1)/iminfo(6)*2 & atsNucAll(:,1) < i/iminfo(6)*2 & atsNucAll(:,8) > 0,9)/divider;
    
    if isempty(cpool)
        cpool = 0;
    end
    mLine(i) = mean(cpool);
    seL(i) = std(cpool)/sqrt(length(cpool));
end

% boundtop = 


figure('pos', [300 200 350 500])
hold off
plot(ATSdistDEmat(:,1)*iminfo(6), ATSdistDEmat(:,2), 'k.', 'markersize', 15);     % L4: x2
axis([1 60  0 y_max ])
xlabel('Radius of Nucleus (\mu\itm\rm)' , 'fontsize',15);
ylabel('DAPI signal intensity (\ita.u.\rm)', 'fontsize',15); 
% text(1.2, 10*1e5, strcat('r = ',num2str(round(corrAll,2,'significant'))), 'color', 'r', 'fontsize', 20);
hold on
plot(1:2:60, mLine, 'c', 'linewidth', 2);
ylabel('Summed ATS intensities (a.u.)' , 'fontsize',15);
xlabel('\itum\rm from distal end', 'fontsize',15);





