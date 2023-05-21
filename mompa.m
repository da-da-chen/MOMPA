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
function [fit,sub_IGD,P_1] =mompa(FUN, varargin)

if ~mod(nargin, 2)
    error('MATLAB:narginchk:notEnoughInputs', ...
        'I have no idea about this, you can guess it');
end
%% Parameter processing
for ind = 1:2:nargin-1
    switch lower(varargin{ind})
        case 'lb'
            lb = varargin{ind + 1};
        case 'ub'
            ub = varargin{ind + 1};
        case 'numobj'
            numObj = varargin{ind + 1};
        case 'max_iter'
            Max_iter = varargin{ind + 1};
        case 'searchagents_no'
            SearchAgents_no = varargin{ind + 1};
        case 'dim'
            dim = varargin{ind + 1};
        case 'numgroup'
            numGroup = varargin{ind + 1};
        case 'minmax'
            if strcmp(varargin{ind + 1}, 'min')
                flagMinMax = 0;
            else
                flagMinMax = 1;
            end
        case 'plotflag'
            plotFlag = varargin{ind + 1};
        otherwise
            error('The function don''t support this parameter');
    end
end
%% Initialization
Iter=0;
FADs=0.2;
P=0.5;
sub_IGD=[];
if flagMinMax
    fit = -inf*ones(SearchAgents_no, numObj);
else
    fit = inf*ones(SearchAgents_no, numObj);
end
%% interative optimization start
if (plotFlag)
    figure(1); hold on;
    h = animatedline;
    h.LineStyle = 'none'; h.Marker = '.'; h.Color = 'r';
    title(FUN)
end

while Iter<Max_iter

    %------------------------------------------------------------
    if Iter==0
        Prey=mompa_initialization(SearchAgents_no,dim,ub,lb);

        [Z,SearchAgents_no] = UniformPoint(SearchAgents_no,numObj);
        [subfit, ~,~] = mompa_getMOFcn(FUN, Prey, numObj);
        Zmin  = min(subfit,[],1);
        Prey=EnvironmentalSelection(FUN,Prey,SearchAgents_no,numObj,Z,Zmin);
        Prey_old=Prey;

    end
    Xmin=repmat(ones(1,dim).*lb,SearchAgents_no,1);
    Xmax=repmat(ones(1,dim).*ub,SearchAgents_no,1);

    Elite = repmat(Prey_old(randi(SearchAgents_no),:),SearchAgents_no,1);

    %------------------------------------------------------------

    CF=(1-Iter/Max_iter)^(2*Iter/Max_iter);
    RL=0.05*mompa_levy(SearchAgents_no,dim,1.5);   %Levy random number vector
    RB=randn(SearchAgents_no,dim);          %Brownian random number vector
    stepsize=zeros(SearchAgents_no,dim);
    for i=1:size(Prey,1)
        for j=1:size(Prey,2)
            R=rand();
            %------------------ Phase 1 (Eq.12) -------------------
            if Iter<Max_iter/3
                stepsize(i,j)=RB(i,j)*(Elite(i,j)-RB(i,j)*Prey(i,j));
                Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j);

                %--------------- Phase 2 (Eqs. 13 & 14)----------------
            elseif Iter>Max_iter/3 && Iter<2*Max_iter/3

                if i>size(Prey,1)/2
                    stepsize(i,j)=RB(i,j)*(RB(i,j)*Elite(i,j)-Prey(i,j));
                    Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j);
                else
                    stepsize(i,j)=RL(i,j)*(Elite(i,j)-RL(i,j)*Prey(i,j));
                    Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j);
                end

                %----------------- Phase 3 (Eq. 15)-------------------
            else

                stepsize(i,j)=RL(i,j)*(RL(i,j)*Elite(i,j)-Prey(i,j));
                Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j);

            end
        end
    end

    %------------------------------------------------------------
    for i=1:size(Prey,1)

        Flag4ub=Prey(i,:)>ub;
        Flag4lb=Prey(i,:)<lb;
        Prey(i,:)=(Prey(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

    end
    Prey_evo=Prey;

    if rand()<FADs
        U=rand(SearchAgents_no,dim)<FADs;
        Prey=Prey+CF*((Xmin+rand(SearchAgents_no,dim).*(Xmax-Xmin)).*U);

    else
        r=rand();  Rs=size(Prey,1);
        stepsize=(FADs*(1-r)+r)*(Prey(randperm(Rs),:)-Prey(randperm(Rs),:));
        Prey=Prey+stepsize;
    end

    for i=1:size(Prey,1)

        Flag4ub=Prey(i,:)>ub;
        Flag4lb=Prey(i,:)<lb;
        Prey(i,:)=(Prey(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

    end
    Prey_fads=Prey;

    Prey_gau = mompa_Gaussian_search(Prey,SearchAgents_no,dim,ones(1,dim).*ub,ones(1,dim).*lb);

    two_Prey = [Prey_old;Prey_evo;Prey_fads;Prey_gau];
    two_loc=unique(two_Prey,'rows');
    [subfit,~ ,P_1] = mompa_getMOFcn(FUN, two_loc, numObj);
    Zmin       = min([Zmin;subfit],[],1);
    Prey=EnvironmentalSelection(FUN,two_loc,SearchAgents_no,numObj,Z,Zmin);
    [fit, ~,~] = mompa_getMOFcn(FUN, Prey, numObj);
    Prey_old=Prey;
    sub_IGD(Iter+1)=mompa_IGD(fit,P_1);

    if (plotFlag && mod(Iter, 1)==0)
        clearpoints(h);
        if numObj==3
            addpoints(h, fit(:, 1), fit(:, 2), fit(:, 3));
        elseif numObj==2
            addpoints(h, fit(:, 1), fit(:, 2));
        end
        drawnow
        %         pause()
    end

    Iter=Iter+1;
    Iter
end






