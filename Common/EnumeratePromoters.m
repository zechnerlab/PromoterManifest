function promoterNames = EnumeratePromoters(folderName)

folderContent = dir(folderName);

counter = 1;
for k=1:length(folderContent)
    if (folderContent(k).isdir == 1)
        
        if (folderContent(k).name(1) ~= '.')
            
            promoterNames{counter} = folderContent(k).name;
            counter = counter + 1;
        end
    end
end