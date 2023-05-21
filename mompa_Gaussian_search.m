% Please refer to the main paper:
% MOMPA: a high performance multi-objective optimizer based on marine predator algorithm
% Long Chen, Fangyi Xu, Kezhong Jin and Zhenzhou Tang
% GECCO '21: Proceedings of the Genetic and Evolutionary Computation Conference Companion
% DOI: https://doi.org/10.1145/3449726.3459581
%        AND
% Marine Predators Algorithm: A nature-inspired metaheuristic
% Afshin Faramarzi, Mohammad Heidarinejad, Seyedali Mirjalili, Amir H. Gandomi
% Expert Systems with Applications
% DOI: https://doi.org/10.1016/j.eswa.2020.113377
% _____________________________________________________
function Prey = mompa_Gaussian_search(Prey,SearchAgents_no,dim,ub,lb)
    for i=1:SearchAgents_no
        d=randperm(dim,1);
        Prey(i,d)=Prey(i,d)+(ub(d)-lb(d))*randn;
    end

    for i=1:size(Prey,1)  
            Flag4ub=Prey(i,:)>ub;
            Flag4lb=Prey(i,:)<lb;    
            Prey(i,:)=(Prey(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
    end 

end