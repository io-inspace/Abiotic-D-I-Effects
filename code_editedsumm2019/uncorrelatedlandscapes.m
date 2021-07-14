%edited 10-18

function [outA, outB, outcorrelationmat,outfractalD_B,outfractalD_A] = uncorrelatedlandscapes(stopper1, stopper2, dodiagnostics,Dcutoff,dimension,iterations,fractaldims,fractalD,landscapes,STD,sa,norm,dim)
saveit = 0;
while stopper1 == 0

%get a landscape B:
    [outlandscapes,outfractalD] = createlandscapes(dimension,iterations,fractaldims,fractalD,landscapes,STD,sa,norm,dodiagnostics,Dcutoff);
    outfractalD_B = outfractalD;

B = cell2mat(outlandscapes);
clear outlandscapes
%get the mean and variance for B
[CmeanB,CvarB] = computemeanvar(dim,B);

%CHECK AND MAKE SURE IT'S NOT CORRELATED WITH ITSELF:
Bcorrmat = zeros(1,length(CmeanB));
for s = 1:length(CmeanB)
 x = reshape(CmeanB{s},1,length(CmeanB{s})^2);
 y = reshape(CvarB{s},1,length(CvarB{s})^2);
 correlation = corrcoef(x,y);
 Bcorrmat(1,s)=correlation(2,1); 
end

if sum(abs(Bcorrmat)>.2)>0 %need a NEW B.
    
else %continue on with A
    
stopper1 = 1;

while stopper2 == 0
   
    [outlandscapes,outfractalD] = createlandscapes(dimension,iterations,fractaldims,fractalD,landscapes,STD,sa,norm,dodiagnostics,Dcutoff);
    outfractalD_A = outfractalD;

    A = cell2mat(outlandscapes);

    [CmeanA,CvarA] = computemeanvar(dim,A);
    [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);

%now that you have A, compare it to itself and B and make sure they are
%uncorrelated, so that any results are not due to a hidden correlation.

    if sum(sum(abs(correlationmat(:,2:3))>.2))>0 %this means that landscapes 
                                             %are correlated with themselves 
                                             %so you need a new A.
    else %if not, go ahead:
        
        if sum(sum(abs(correlationmat) > .2))==0 %if all correlations are under .2, stop the loop you're done.
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else 
             Arot1 = rot90(A);
             [CmeanA,CvarA] = computemeanvar(dim,Arot1);
             [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Arot1;
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else
                Arot2 = rot90(Arot1);
                [CmeanA,CvarA] = computemeanvar(dim,Arot2);
                [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Arot2;
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else        
                Arot3 = rot90(Arot2);
                [CmeanA,CvarA] = computemeanvar(dim,Arot3);
                [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Arot3; 
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else
                Arot4 = rot90(Arot3);
                [CmeanA,CvarA] = computemeanvar(dim,Arot4);
                [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Arot4;
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else
                Aud1 = rot90(Arot1);
                [CmeanA,CvarA] = computemeanvar(dim,Aud1);
                [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Aud1;
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else
                Aud2 = rot90(Arot2);
                [CmeanA,CvarA] = computemeanvar(dim,Aud2);
                [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Aud2;
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else
                Aud3 = rot90(Arot3);
                [CmeanA,CvarA] = computemeanvar(dim,Aud3);
                [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Aud3;
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else
                Aud4 = rot90(Arot4);
                [CmeanA,CvarA] = computemeanvar(dim,Aud4);
                [correlationmat] = corrsfunction(CmeanA,CmeanB,CvarA,CvarB);
        end
        if sum(sum(abs(correlationmat) > .2))==0
            A = Aud4;
            stopper2 = 1;
            if saveit==1
                filename = ['uncorrelatedlandscapes_dim' num2str(dim) '.mat'];
                save(filename,'A','B','correlationmat')
            end
        break
        else
        stopper2 = 0;
        end
    end
end

end

end

outA = A;
outB = B;
outcorrelationmat = correlationmat;

end
