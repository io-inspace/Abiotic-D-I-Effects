%before, when I did the oldind matrix, all newly germinated perennials 
%would survive, which is stupid.  They're subject to mortality.  I don't
%know whether I want them to produce seeds or not.

function [outputnewspsmatrix] = mortality(dim,spsmatrix,maxsp,perennials,mort,gen)
rng(gen*0.6803); %because I'm going to use gen to seed this in multiple fxns, I multiplied by a randomly generated number.
newspsmatrix = zeros(dim,dim,maxsp);
 %annuals: they all die, so it'll be all 0s anyway
 %perennials:
   
   %next, a loop through each species
    for sps = 1:length(perennials)
        occupied = find(spsmatrix(:,:,perennials(sps))>0);
        
            for o = 1:length(occupied); %do this microsite by microsite
                %get the x,y coordinate for the occupied(o) you're on.
                getdim = size(spsmatrix);
                [x,y]=ind2sub(getdim(1:2),occupied(o));
                
                numold = spsmatrix(x,y,perennials(sps)); %how many individuals are there
                liveordie = rand(1,numold); %draw a random number for each of these to see if they live or die
                    for decide = 1:length(liveordie) %go through each of these individuals
                        if liveordie(decide) > mort %if they draw a lucky number, they survive to the next gen
                            newspsmatrix(x,y,perennials(sps)) = newspsmatrix(x,y,perennials(sps))+1;
                        end
                    end
            end
    end

outputnewspsmatrix = newspsmatrix;

end
