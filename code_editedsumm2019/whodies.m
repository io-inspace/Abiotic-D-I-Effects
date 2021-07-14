
function [newentry] = whodies(whichsps,normedprob,germind_entry,remove,spsmatrix_entry,openspace,k)

helpfulmatrix = [whichsps,normedprob']; %get sps and their probability of dying
            if any(helpfulmatrix(:,2) > .9999)
              %then remove all of the species that have no chance of surviving first
              %and then whoever is left.
              nogerm = helpfulmatrix(helpfulmatrix(:,2)>.9999,1);%this is the sps that for sure dies.
              if germind_entry(:,:,nogerm) <= remove 
              germind_entry(:,:,nogerm) = 0;
              else
                  germind_entry(:,:,nogerm) = germind_entry(:,:,nogerm)-remove;
              end
              newentry = spsmatrix_entry+germind_entry;
              if sum(newentry) > k
                  yesgerm = helpfulmatrix(helpfulmatrix(:,2)<.9999,1);
                  newentry(:,:,yesgerm) = k;
              end
            else
                    previous = helpfulmatrix(1,2);
                    %loop of death
                    %remove = sum(germind(x,y,:))-openspace; %how many individuals need to be removed
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
                       germind_entry(:,:,finalprobs(1,1))=germind_entry(:,:,finalprobs(1,1))-1; 
                     else 
                     for w = 1:(length(column)-1)
                       if column(w) < fate && fate < column(w+1)
                         germind_entry(:,:,finalprobs(w+1,1))=germind_entry(:,:,finalprobs(w+1,1))-1;
                       end
                     end
                     end
                     newentry = spsmatrix_entry+germind_entry;
                    else %if remove > 1 
                      %densities = vertcat(germind(x,y,whichsps));
                      densities = vertcat(germind_entry(:,:,whichsps));
                      finalprobs = horzcat(helpfulmatrix(:,1),column,densities(:));
                      while  sum(germind_entry(:,:,:))>openspace %how many rounds of "who doesn't make it" are we playing?
                      fate = rand(1,1);
                            if fate < finalprobs(1,2) && finalprobs(1,3)>0;
                              germind_entry(:,:,finalprobs(1,1))=germind_entry(:,:,finalprobs(1,1))-1; 
                              densities = vertcat(germind_entry(:,:,whichsps));
                              finalprobs = horzcat(helpfulmatrix(:,1),column,densities(:));
                            else  
                                for w = 1:(length(column)-1)
                                  if column(w) < fate && fate < column(w+1) && finalprobs(w+1,3)>0
                                     germind_entry(:,:,finalprobs(w+1,1))=germind_entry(:,:,finalprobs(w+1,1))-1;
                                     densities = vertcat(germind_entry(:,:,whichsps));
                                     finalprobs = horzcat(helpfulmatrix(:,1),column,densities(:));
                                  end
                                end
                            end
                      end
                      newentry = spsmatrix_entry+germind_entry;
                    end
           end
end