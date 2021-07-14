for ls = 1:20
saveit = 1;

load(['dim_65_' num2str(ls) '_neutRdiffG#diffC_uni_1_5.mat'])
jitter = .9;


all_locations_bysps = cell(4,1);
for sp = 1:4
total=endallinds{1,1}(:,:,sp);

allinds=sum(sum(total,2));

occupied = find(total>0);
all_locations = cell(length(occupied),sp); 
for o = 1:length(occupied); %do this microsite by microsite
    %get the x,y coordinate for the occupied(o) you're on.
    getdim = size(total);
    [x,y]=ind2sub(getdim(1:2),occupied(o));
    %how many individuals are at this location?
    numinds = zeros(total(x,y),3);
            j1 = -(jitter);
            j2 = jitter;%(you have to do +1 with floor)
            jitterit = (j2-j1).*rand(total(x,y),2)+j1; %draw a random location
            
        for j = 1:total(x,y)
            numinds(j,1)=sp;
            numinds(j,2)= (y+jitterit(j,2)); %this is so fucked up but it gives the correct map in R.
            numinds(j,3)=65-(x+jitterit(j,1));
        end
        all_locations{o,sp} = numinds;
end
all_locations_bysps{sp,1} = vertcat(all_locations{:});
end

spsmap = vertcat(all_locations_bysps{:});

if(saveit==1)
    spsmap = array2table(spsmap,'VariableNames',{'ID','X','Y'});

    writetable(spsmap,['ppmap_ls' num2str(ls) '.csv'],'Delimiter',',','QuoteStrings',true)

end
end
