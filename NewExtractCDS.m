function [CDS] = NewExtractCDS(gbfile,Species)


%%
switch Species 
    
    case 1 %S. cerevisiae 
        GenBankFile = genbankread(gbfile);
        CDS = featureparse(GenBankFile,'feature','cds','Sequence','true');
        
        %Just for S. cerevisiae - I need the dbx_ref to link to the SGD for when I group genes according to ontology
        for i=1:height(struct2table(CDS))
        CDS(i).db_xref = CDS(i).db_xref{1,2};
        end 
        
        
        
    case 2 %S. pombe
        
        GenBankFile = genbankread(gbfile);
        CDS = featureparse(GenBankFile,'feature','cds','Sequence','true');
       
        
    case 3 %A. gossypii
    
        GenBankFile = genbankread(gbfile);
        CDS = featureparse(GenBankFile,'feature','cds','Sequence','true');
        
    
end 

        CDS = struct2table(CDS);
        
        for i=height(CDS):-1:1
            
            Sequence = string(extractBetween(CDS.Sequence(i), 1, strlength(string(CDS.Sequence(i)))-3)); 
            tag_number = strfind(Sequence,'tag');
            taa_number = strfind(Sequence,'taa');
            tga_number = strfind(Sequence,'tga'); 
            stopcodon_occurrences = [tag_number taa_number tga_number] - 1;
            idx = find(mod(stopcodon_occurrences,3)==0);
            
            if ~isempty(idx) || strlength(CDS.Sequence(i)) < 300 || strcmp(CDS.codon_start(i),'1') == 0 || rem(strlength(string(CDS.Sequence(i))),3) ~= 0
                CDS(i,:) = [];
            end  
            
            
        end
        
        Strand = strings(height(CDS),1);
        
        for i=1:height(CDS)
            if extractBetween(string(CDS.Location(i)),1,4) == 'comp'
                Strand(i,1) = 'reverse';
            else Strand(i,1) = 'forward';
            end 
        end 
        CDS.Strand = Strand; 
        CDS = CDS(:,{'locus_tag', 'protein_id', 'db_xref', 'translation', 'Sequence', 'Strand'});

end

