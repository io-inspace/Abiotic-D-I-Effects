%edited 10-18

%diagonals
function [outlandscape] = diagonals(landscape,delta,f,k,sa)
 %MAIN    
    main = diag(landscape,0); %pull out the main diagonal and call it main
    points = find(main); %find the non-zeros: these are the numbers we need to calculate the midpoints (mps)
  for y = 1:length(points)-1; %determine how many mps are needed
    mp = (points(y)+points(y+1))/2; %mp is a location, not a number
    davg = (landscape(mp-f,mp-f)+landscape(mp-f,mp+f)+landscape(mp+f,mp-f)+landscape(mp+f,mp+f))/4; %davg are numbers from the 
                                                                                                    %landscape matrix
    main(mp) = davg + delta*randn; %this is the final entry: the average of the surrounding points + some wiggle
  end
matrix = diag(main,0); %stick this back into a temporary matrix called matrix, specifically for the diagonals of landscape.

%UPPER TRIANGLE 
    
  for x = 1:k; %go diagonal by diagonal
    if any(diag(landscape,x)) == 1; %check to see if there are nonzero elements in the diagonal.
       diagonal = diag(landscape,x); %if there are, pull it out
       points = find(diagonal); %now follow the same process as for the main diagonal
    for y = 1:length(points)-1;
       mp = (points(y)+points(y+1))/2;
       ump = mp + x; %this is to account for the fact that we're in the upper triangle, so the distance to neighbors f needs 
                     %to be modified by an amount mp + x.  For example, in a 5x5 matrix landscape(mp-f,ump-f) is (1,3), which
                     %puts us in the upper triangle for the top left point of the box, whereas without ump we'd get (1,1)
       davg = (landscape(mp-f,ump-f)+landscape(mp+f,ump-f)+landscape(mp-f,ump+f)+landscape(mp+f,ump+f))/4;
       diagonal(mp) = davg + delta*randn;
    end 
    matrix = matrix + diag(diagonal,x); %put the diagonal into matrix with the main
    else %skip those diagonals that are currently empty
        diagonal = diag(landscape,x);
    matrix = matrix + diag(diagonal,x);
    end
  end
  
%LOWER TRIANGLE
%repeat the same process as the upper triangle
    
  for x = 1:k;
    if any(diag(landscape,-x)) == 1;
       diagonal = diag(landscape,-x);
       points = find(diagonal);
    for y = 1:length(points)-1;
       mp = (points(y)+points(y+1))/2;
       lmp = mp + x; %equivalent to ump but for the lower triangle.
       davg = (landscape(lmp-f,mp-f)+landscape(lmp+f,mp-f)+landscape(lmp+f,mp+f)+landscape(lmp-f,mp+f))/4;
       diagonal(mp) = davg + delta*randn;
    end 
    matrix = matrix + diag(diagonal,-x);
    else
        diagonal = diag(landscape,-x);
    matrix = matrix + diag(diagonal,-x);
    end
  end

landscape = matrix; %rename the newly constructed matrix of diagonals and proceed.

%do random additions:
if sa == 1;
for i = 1:length(find(landscape));%this goes through each non-zero entry
    totpnts = find(landscape);%pull out the non-zeros
    landscape(totpnts(i)) = landscape(totpnts(i))+delta*randn;%and add some wiggle to each of them
end
end

outlandscape = landscape;

end
