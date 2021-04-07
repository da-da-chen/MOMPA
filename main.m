warning('off')
clear, clc
a_mompa=[];
 for i=1:1
%%
    %ZDT系列
        %           M    D	Range                   
        %ZDT1-ZDT3	2	30	xi∈[0,1],1≤i≤D	
        %ZDT4       2	10	x1∈[0,1],xi∈[-5,5],2≤i≤D	'lb', [0 -5*ones(1,9)], 'ub', [1 5*ones(1,9)]
        %ZDT6       2	10	xi∈[0,1],1≤i≤D	
    %DTLZ系列
        %              M   D	Range  
        %DTLZ1          3	7	xi∈[0,1],1≤i≤D	
        %DTLZ2-DTLZ6	3	12	xi∈[0,1],1≤i≤D	
        %DTLZ7          3	22	xi∈[0,1],1≤i≤D	

%     [fit,IGD,P] =   mompa('DTLZ4', 'Max_iter', 3000, 'SearchAgents_no',100,...
%             'minmax', 'min', 'plotFlag', 1, 'dim', 12, 'numObj', 3, ...
%             'numgroup', 1, 'lb', 0, 'ub', 1);
%         i
%         a_mompa(i) = IGD(end);

%%
    %WFG系列
    %          M   D	Range  
    %WFG3-WFG9 3  12	zi∈[0,2i],1≤i≤D	

    final_lb =zeros(1,12);
    final_ub =2 : 2 : 2*12;%
    [fit,IGD,P] = mompa('WFG9', 'Max_iter', 3000, 'SearchAgents_no',100,...
            'minmax', 'min', 'plotFlag', 1, 'dim', 12, 'numObj', 3, ...
            'numgroup', 1, 'lb', final_lb, 'ub', final_ub);
         i
         a_mompa(i) = IGD(end);
 end
a_mompa
b1=mean(a_mompa);
b2=std(a_mompa);


scatter3(P(:, 1), P(:, 2),P(:,3));
 %scatter(P(:, 1), P(:, 2));
% grid on;  