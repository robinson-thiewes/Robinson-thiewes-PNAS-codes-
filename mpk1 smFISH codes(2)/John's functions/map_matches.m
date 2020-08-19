function matched_rna = map_matches(prna,rnas_atsEx, ExATSrna_ATS, ChOrder, iminfo)
%   output: matches -- struct of RNAs that were mapped to each other.
%   Contains 8 fields:
%   ID: RNA ID (example intron 46 would be "in46", exon 12 would be "ex12")
%   | location_um: RNA location in micrometers [x,y,z] | match_ID: matched 
%   RNA ID | match_location_um: matched RNA location in micrometers 
%   [x, y, z] | distance_um: distance between matched RNAs in micrometers |
%   match_intensity: matched RNA intensity | location: RNA location in 
%   pixels [x,y,z] | match_location: matched RNA location in pixels [x,y,z]


    temp_intron = prna{ChOrder == 1}((prna{ChOrder == 1}(:,9) ~= -1),:);
    temp_exon = ExATSrna_ATS((ExATSrna_ATS(:,9) ~= -1),:);
    
    matched_rna = struct('ID', {}, 'location_um', {}, 'match_ID', ...
        {}, 'match_location_um', {},'distance_um',{}, 'match_intensity', {},...
        'location', {}, 'match_location', {});
    for i = [1:size(temp_intron,1)]
        matched_rna(i).ID = strcat("in",string(temp_intron(i,8)));
        matched_rna(i).location = temp_intron(i,[1,2,7]);
        
        matched_rna(i).match_ID = strcat("ex",string(temp_intron(i,9)));
        matched_rna(i).match_location = rnas_atsEx((rnas_atsEx(:,8) == temp_intron(i,9)),[1,2,7]);
        matched_rna(i).match_intensity = temp_intron(i,10);
    end
    exon_num = 0;
    for i = [size(temp_intron,1)+1:size(temp_intron,1)+size(temp_exon,1)]
        exon_num = exon_num + 1;
        %fprintf('i: %d | exon_num: %f\n',i,exon_num)
        matched_rna(i).ID = strcat("ex",string(temp_exon(exon_num,8)));
        matched_rna(i).location = temp_exon(exon_num,[1,2,7]);
        matched_rna(i).match_ID = strcat("in",string(temp_exon(exon_num,9)));
        matched_rna(i).match_location = prna{ChOrder == 1}...
            ((prna{ChOrder == 1}(:,8) == temp_exon(exon_num,9)),[1,2,7]);
        matched_rna(i).match_intensity = temp_exon(exon_num,10);
    end
    
    for i = 1:size(matched_rna,2)
        matched_rna(i).location_um = matched_rna(i).location .* iminfo(5:7).';
        matched_rna(i).match_location_um = matched_rna(i).match_location .* iminfo(5:7).';
        distance_difference_um = matched_rna(i).match_location_um - matched_rna(i).location_um;
        matched_rna(i).distance_um = sqrt(sum(distance_difference_um .^2));
    end
end