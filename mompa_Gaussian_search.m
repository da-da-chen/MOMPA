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