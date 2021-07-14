%returns species richness (should I also do variation in richness?) at each
%scale.
function [outaverageSoverIC] = averagerichness(endrichness,maxsp,scale,initcond)

presabsdata = cell(maxsp,length(scale),initcond);
abundancedata = endrichness;
        for initc = 1:initcond
            for sampunit = 1:length(scale)
                for species = 1:maxsp
                    presabsdata{species,sampunit,initc} = zeros(size(abundancedata{species,sampunit,initc}));
                    where = find(abundancedata{species,sampunit,initc});
                    for quadrat = 1:length(where);
                        if abundancedata{species,sampunit,initc}(where(quadrat)) > 1
                            presabsdata{species,sampunit,initc}(where(quadrat)) = 1;
                        end
                    end
                end
            end
        end
        totalS = cell(1,length(scale),initcond);
        for first = 1:initcond
            for second = 1:length(scale)
                richindex = presabsdata{1,second,first};
                for third = 2:maxsp
                    richindex = plus(richindex,presabsdata{third,second,first});
                end
                totalS{1,second,first} = richindex;
            end
        end
        averageSoverIC = cell(1,length(scale));
        for second = 1:length(scale)
            totalSindex = totalS{1,second,1};
            for first = 2:initcond
                totalSindex = totalSindex+totalS{1,second,first};
            end
            averageSoverIC{1,second} = totalSindex/initcond;
        end
        
outaverageSoverIC = averageSoverIC;
        
end
