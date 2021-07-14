%EDITED 1-26

function [outputmatrix] = competition_sym(spsmatrix,total,k,alpha,sigsq,C,npC,diffC,savg_comps)
overk = find(total>k);%gives you the indices to these microsites that have occupancy > k
rng(2);
for comp = 1:length(overk)%do this microsite by microsite
    
    %get the x,y coordinate for the overk(comp) you're on.
    getdim = size(spsmatrix);
    [x,y]=ind2sub(getdim(1:2),overk(comp));
    
    whichsps = find(spsmatrix(x,y,:)); %which sps are in the microsite
            
    %find the R(x,y) of each species - this is the competitive
    %ability of each species defined by a gaussian fxn.  The closer
    %a species niche (npC) is to C(x,y), the higher its R.
    
    %ugh now this is all messed up because you don't want the second one in
    %the diffC sims.
    Rbysps = zeros(1,length(whichsps));
        for bysps = 1:length(whichsps);
            if diffC == 1 %with spatial variation in competitive ability
            Rbysps(1,bysps) = alpha*exp((-.5)*(C(x,y)-(npC(whichsps(bysps))))^2/(sigsq));
            else %without spatial variation in competitive ability
            Rbysps(1,bysps) = savg_comps(1,whichsps(bysps));
            end
        end
       
        
            
    %find r of species j, which is the probability it dies due to
    %competition.
    probabilities = zeros(1,length(whichsps));
       for j = 1:length(whichsps); 
          compeffect = zeros(1,length(whichsps));
          for m = 1:length(whichsps); 
             if j == m
               compeffect(m) = 0;
             else  
                 compeffect(m) = Rbysps(m)*alpha*spsmatrix(x,y,whichsps(m));%and then alpha again
             end
          end
           probabilities(j) = spsmatrix(x,y,whichsps(j))*(spsmatrix(x,y,whichsps(j))+sum(compeffect));
       end
        %normalize so that all values sum to 1:
        normedprob = probabilities/sum(probabilities);
           
   %randomly draw who dies based on these probabilities.          
   helpfulmatrix = [whichsps,normedprob'];
   %loop of death
   remove = total(overk(comp)) - k; 
   previous = helpfulmatrix(1,2);
   column = zeros(1,(size(helpfulmatrix,1)-1));
   for stack = 1:(size(helpfulmatrix,1)-1)
     new = previous+helpfulmatrix(1+stack,2);
     column(1,stack) = new;
     previous = new;
   end
     column = [helpfulmatrix(1,2),column]';
     finalprobs = horzcat(helpfulmatrix(:,1),column); 
   if remove == 1
     fate = rand(1,1);
     if fate < finalprobs(1,2);
       spsmatrix(x,y,finalprobs(1,1))=spsmatrix(x,y,finalprobs(1,1))-1; 
     else 
     for w = 1:(length(column)-1)
       if column(w) < fate && fate < column(w+1)
         spsmatrix(x,y,finalprobs(w+1,1))=spsmatrix(x,y,finalprobs(w+1,1))-1;
       end
     end
     end
   else %if remove > 1 
      densities = vertcat(spsmatrix(x,y,whichsps));
      finalprobs = horzcat(helpfulmatrix(:,1),column,densities(:));
      while  sum(spsmatrix(x,y,:))>k %how many rounds of "who dies" are we playing?
      fate = rand(1,1);
            if fate < finalprobs(1,2) && finalprobs(1,3)>0;
              spsmatrix(x,y,finalprobs(1,1))=spsmatrix(x,y,finalprobs(1,1))-1; 
              densities = vertcat(spsmatrix(x,y,whichsps));
              finalprobs = horzcat(helpfulmatrix(:,1),column,densities(:));
            else  
                for w = 1:(length(column)-1)
                  if column(w) < fate && fate < column(w+1) && finalprobs(w+1,3)>0
                     spsmatrix(x,y,finalprobs(w+1,1))=spsmatrix(x,y,finalprobs(w+1,1))-1;
                     densities = vertcat(spsmatrix(x,y,whichsps));
                     finalprobs = horzcat(helpfulmatrix(:,1),column,densities(:));
                  end
                end
            end
      end
   end
end
outputmatrix = spsmatrix;
end
