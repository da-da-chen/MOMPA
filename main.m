% Developed in MATLAB R2019b
% Source codes demo version 1.0
% _____________________________________________________
%  Author, inventor and programmer of MOMPA: Long Chen,
%  Email: chenlong@zjnu.edu.cn
%  Co-authors:Xuebing Cai, Kezhong Jin, Zhenzhou Tang
% _____________________________________________________
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off')
clear, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      ZDT
%Parameters      numObj      dim	     lb & ub
%ZDT1-ZDT3	       2	     30	    xi∈[0,1],1≤i≤dim
%                                   final_lb = 0;
%                                   final_ub = 1;
%ZDT4              2	     10	    x1∈[0,1],xi∈[-5,5],2≤i≤dim
%                                   final_lb = [0 -5*ones(1,9)];
%                                   final_ub = [1 5*ones(1,9)];
%ZDT6              2	     10	    xi∈[0,1],1≤i≤dim
%                                   final_lb = 0;
%                                   final_ub = 1;
%%                      DTLZ
%Parameters      numObj      dim	     lb & ub
%DTLZ1             3	     7	    xi∈[0,1],1≤i≤dim
%                                   final_lb = 0;
%                                   final_ub = 1;
%DTLZ2-DTLZ6	   3	     12	    xi∈[0,1],1≤i≤dim
%                                   final_lb = 0;
%                                   final_ub = 1;
%DTLZ7             3	     22	    xi∈[0,1],1≤i≤dim
%                                   final_lb = 0;
%                                   final_ub = 1;
%%                       WFG
%Parameters      numObj      dim	     lb & ub
%WFG3-WFG9         3         12	    xi∈[0,2i],1≤i≤dim
%                                   final_lb = zeros(1,12);
%                                   final_ub = 2 : 2 : 2*12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_mompa=[];
for i=1:1
        final_lb =0;      %Different test functions choose different parameters
        final_ub =1;      %Different test functions choose different parameters
        [fit,IGD,P] = mompa('DTLZ1', ...  %test function name
        'Max_iter', 3000, ...  
        'SearchAgents_no',100,... 
        'minmax', 'min', ...      
        'plotFlag', 1, ...        
        'dim', 7, ...      %  Different test functions choose different parameters      
        'numObj', 3, ...    %  Different test functions choose different parameters  
        'numgroup', 1, ...        
        'lb', final_lb, ...       
        'ub', final_ub);         
    i
    a_mompa(i) = IGD(end);
end
a_mompa
b1=mean(a_mompa);
b2=std(a_mompa);

if size(P,2)==3
    scatter3(P(:, 1), P(:, 2),P(:,3));
else 
    scatter(P(:, 1), P(:, 2));
end