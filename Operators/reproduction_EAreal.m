function QWE=reproduction_EAreal(Pop, EP, D, lu)
    % reproduction operation to produce children
    
    % mating selection is needed for most operators
    if rand <0.5
        ind1=EP(randi([1, size(EP,1)]),:); % choose one from Archive
        while true
            ind2_id= tournament_selection(2,1, Pop(:,end-1));% choose one from current population
            ind2=Pop(ind2_id,1:end-1);
            
            if ~all(ind1==ind2), break;end % make sure two selected inds are different
        end
        mating=[ind1; ind2];
    else
        while true
            mating_id=tournament_selection(2,2, Pop(:,end-1));
            mating =Pop(mating_id,1:end-1);

            if ~all(mating(1,:)==mating(2,:)), break;end % make sure two selected inds are different
        end
    end
    

    QWE = EAreal(mating,D,lu); % this generates two children
    QWE = QWE(randi([1,2]),:);
end


function index = tournament_selection(K,N,F)
%TournamentSelection - Tournament selection
%
%   P = TournamentSelection(K,N,fitness) returns the indices
%   of N solutions by K-tournament selection based on their fitness values.
%   In each selection, the candidate having the minimum fitness1 value will
%   be selected.
%
%   Example:
%       P = TournamentSelection(2,100,FrontNo)
    [~,rank] = sortrows(F);
    [~,rank] = sort(rank);
    Parents  = randi(length(F),K,N);
    [~,best] = min(rank(Parents),[],1);
    index    = Parents(best+(0:N-1)*K);
end