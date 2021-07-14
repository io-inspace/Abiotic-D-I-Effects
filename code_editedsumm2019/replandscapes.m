%edited 12-2018

%get landscapes randomly populated with a certain number
%of inidividuals of each species.

%replicating across initial conditions:

function [originallandscapes] = replandscapes(initcond,dim,maxsp)

%define how many individuals total to start with 
Ntot = dim^2; 
%split up the total evenly among all species.
NI = repmat(floor(Ntot/maxsp),1,maxsp); 

originallandscapes = cell(1,initcond);
rng(1);
for ic = 1:initcond
    
%STEP 1: RANDOM PLACEMENT
spsmatrix = zeros(dim,dim,maxsp);

for n = 1:length(NI); %do for each species
    for i = 1:NI(n); %do for each individual
        a = 1;
        b = dim+1;%(you have to do +1 with floor)
        r = floor((b-a).*rand(1,2)+a); %draw a random location
        spsmatrix(r(1),r(2),n) = spsmatrix(r(1),r(2),n)+1;
     end
end

originallandscapes{1,ic} = spsmatrix;

end
