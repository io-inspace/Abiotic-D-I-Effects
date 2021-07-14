%edited 12-11

function [endallinds,endtotalinds,endrichness,endfitness,endprobdeath] = simulationv2018(dist_a,dist_p,R_a,R_p,C,G,dim,npR,npC,npG,generations,initcond,scale,maxsp,sigsq,sigsqR,Ma,Mp,alpha,beta,m,k,prop,perennials,originallandscapes,diffC,savg_comps,savg_seeds,savg_germs,diffG,simR)
%rng(2); %if you try this way, make sure you go through and annotate all the
%other rngs.

%PARAMETERS

%constant mortality
mort = 1/m; 

%OUTPUT
richness = cell(maxsp,length(scale),initcond);
allinds = cell(1,initcond);
fitness = cell(initcond,length(scale));
%inter = cell(initcond,1);
%intra = cell(initcond,1);
probdeath = cell(initcond,1);
totalinds = zeros(generations,maxsp,initcond);

for ic = 1:initcond
    spsmatrix = originallandscapes{1,ic};
    
%STEP 2: INITIAL (SYMMETRIC) COMPETITION

%identify microsites where there are more than k individuals
total = sum(spsmatrix,3); %now you have a matrix telling you how many 
%individuals are in each microsite

%call the function competition_sym that will remove species from microsites with 
%>k individuals based on the probability function.
[outputmatrix]=competition_sym(spsmatrix,total,k,alpha,sigsq,C,npC,diffC,savg_comps);

%this is the spsmatrix used going into the generation loop.
spsmatrix = outputmatrix;

%STEP 3: the big loop (germination, competition, reproduction, and dispersal) for gen generations

for gen = 1:generations
 if gen == 1
          
 %CENSUS INDIVIDUALS
 %total no. individuals on the landscape (by sps):
 for species = 1:maxsp
    totalinds(gen,species,ic) = sum(sum(spsmatrix(:,:,species))); %#annuals+perennials before mortality
 end
 %allinds = cell(generations,initcond);
 %allinds{gen,ic} = spsmatrix;
 total = sum(spsmatrix,3); %how many inidividuals are in a microsite after competition
 
 %REPRODUCTION & DISPERSAL
    [outputseedmatrix] = reproduction_dispersal(dist_a,dist_p,dim,maxsp,total,spsmatrix,Ma,Mp,R_a,R_p,npR,prop,sigsqR,gen,simR,savg_seeds);
    seedind = outputseedmatrix;
    %allseeds = cell(generations,initcond);
    %allseeds{gen,ic} = seedind; 
 else %gen~=1
    if gen == generations
        surviving_perennials = spsmatrix; 
%NOPE, needs to be before dispersal. %this is the "seed input" at the start of the
                                           %final generation from the
                                           %generation before. It includes
                                           %surviving perennials and all
                                           %seeds AFTER DISPERSAL.
        %individuals_in = allinds{generations-1,ic};
        %seedsin = seedoutput(maxsp,dim,individuals_in,Ma,Mp,R,npR,sigsq);
        seedind_IN = surviving_perennials+seedind;
    end
    
     %GERMINATION
     [outputgermind] = germination(dim,maxsp,seedind,G,npG,beta,sigsq,gen,diffG,savg_germs);
     germind = outputgermind;
     
     %Now that we have a bunch of germinated seedlings, record the
     %[[competitive environment]] experienced by each sps.
     if gen == generations
       [dyingprob] = compenv(germind,spsmatrix,dim,maxsp,C,npC,alpha,sigsq,diffC,savg_comps,k);
     end
     %ASYMMETRIC COMPETITION
     [outputcompmatrix] = competition_asym(germind,spsmatrix,dim,k,alpha,sigsq,C,npC,gen,diffC,savg_comps);
     spsmatrix = outputcompmatrix;

     total = sum(spsmatrix,3);

%CENSUS INDIVIDUALS
%total no. individuals on the landscape (by sps):
for species = 1:maxsp
    totalinds(gen,species,ic) = sum(sum(spsmatrix(:,:,species))); %#annuals+perennials before mortality
end



if gen == generations 
    
   allinds{1,ic} = spsmatrix; %individuals before reproduction and dispersal, 
                             %which is what you need to calculate seeds.
                             %This at (generations-1) +
                             %surviving_perennials gives seedind_IN.
   %split this up into two output codes.  First, richness:
   [outputrichness] = calculateS(richness,scale,spsmatrix,dim,maxsp,ic);
   richness = outputrichness;

   %then fitness:
   %at this point, this is just the total seeds output before dispersal in the
   %last generation:
   %[outputfitness] = calculatefit(scale,fitness,allinds,fitls,R,npR,Ma,Mp,dim,maxsp,ic,sigsq); 
   %fitness = outputfitness;
end

%now we WANT THE SEED MATRIX FROM THE FINAL GENERATION
%if gen < generations
    [outputseedmatrix] = reproduction_dispersal(dist_a,dist_p,dim,maxsp,total,spsmatrix,Ma,Mp,R_a,R_p,npR,prop,sigsqR,gen,simR,savg_seeds);
    seedind = outputseedmatrix;
%end

 end
 
%...INCLUDING the surviving perennials.
%if gen < generations 
 %MORTALITY
 
 [outputnewspsmatrix] = mortality(dim,spsmatrix,maxsp,perennials,mort,gen);
 spsmatrix = outputnewspsmatrix; %this has surviving perennials (which should be treated like seeds with new fitness measurements)
%end
if gen == generations
    individuals_in = allinds{1,ic}; %individuals surviving competition before repro+disp
    seedsout_final = seedoutput(maxsp,dim,individuals_in,Ma,Mp,R_a,R_p,npR,sigsqR,gen,simR,savg_seeds);
    seedind_OUT = seedsout_final+spsmatrix;
    seedfrac = seedind_OUT./seedind_IN;
    [fitness] = calculatefit(scale,fitness,seedfrac,dim,maxsp,ic);
end
    
end
%inter{ic,1} = inter_compeffect;
%intra{ic,1} = intra_compeffect;
probdeath{ic,1} = dyingprob; %from compenv.m
end

endrichness = richness;
endallinds = allinds;
endfitness = fitness;
endtotalinds = totalinds;
%endcomp_inter = inter;
%endcomp_intra = intra;
endprobdeath = probdeath;

end

