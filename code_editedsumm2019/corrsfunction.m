%edited 10-18

%it's worth it to write a function for the flipping landscapes around part
function correlationmat = corrsfunction(meanA,meanB,varA,varB)

correlationmat = zeros(length(meanA),6);

for m = 1:length(meanA)
    x = reshape(meanA{m},1,length(meanA{m})^2);
    y = reshape(meanB{m},1,length(meanB{m})^2);
    correlation = corrcoef(x,y);
    correlationmat(m,1)=correlation(2,1); 
    
    x = reshape(meanA{m},1,length(meanA{m})^2);
    y = reshape(varA{m},1,length(varA{m})^2);
    correlation = corrcoef(x,y);
    correlationmat(m,2)=correlation(2,1); 
    
    x = reshape(meanB{m},1,length(meanB{m})^2);
    y = reshape(varB{m},1,length(varB{m})^2);
    correlation = corrcoef(x,y);
    correlationmat(m,3)=correlation(2,1); 
    
    x = reshape(meanA{m},1,length(meanA{m})^2);
    y = reshape(varB{m},1,length(varB{m})^2);
    correlation = corrcoef(x,y);
    correlationmat(m,4)=correlation(2,1); 
    
    x = reshape(meanB{m},1,length(meanB{m})^2);
    y = reshape(varA{m},1,length(varA{m})^2);
    correlation = corrcoef(x,y);
    correlationmat(m,5)=correlation(2,1); 
    
    x = reshape(varA{m},1,length(varA{m})^2);
    y = reshape(varB{m},1,length(varB{m})^2);
    correlation = corrcoef(x,y);
    correlationmat(m,6)=correlation(2,1); 
end   
end