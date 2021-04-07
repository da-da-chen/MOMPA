function Population= EnvironmentalSelection(FUN,Population,N,M,Z,Zmin)
% The environmental selection of NSGA-III

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end

    %% Non-dominated sorting
    [Population_objs,~ ]= mompa_getMOFcn(FUN, Population, M);
    [FrontNo,MaxFNo] = NDSort(Population_objs,N);
    Next = FrontNo < MaxFNo;%进入下一代
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);%==最后一层的下标
    Choose = LastSelection(Population_objs(Next,:),Population_objs(Last,:),N-sum(Next),Z,Zmin);%基于参考点后最后一层选择k个 逻辑值0、1
    Next(Last(Choose)) = true;%
    % Population for next generation
    Population = Population(Next,:);
%     Population_objs = getMOFcn(FUN, Population, M);
%     [FrontNo,~] = NDSort(Population_objs,N);
%     sort_No_y=sortrows([FrontNo',Population_objs ,Population],[1,2]);
    
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);%减去理想点
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);%A的逆矩阵是 A1，则 B/A = B*A1，A\B = A1*B
    a = 1./Hyperplane;
    if any(isnan(a))%若a的元素为NaN（非数值），在对应位置上返回逻辑1（真），否则返回逻辑0（假）
        a = max(PopObj,[],1)';
    end
    % Normalization
    %a = a-Zmin';
    PopObj = PopObj./repmat(a',N,1);%如果a、b是矩阵，a./b就是a、b中对应的每个元素相除，得到一个新的矩阵；

    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');%cosine’ 针对向量,夹角余弦距离Cosine distance(‘cosine’)余弦相似度= 1-余弦距离

    %杰卡德距离Jaccard distance(‘jaccard’)Jaccard距离常用来处理仅包含非对称的二元(0-1)属性的对象。很显然，Jaccard距离不关心0-0匹配[1]。
    %夹角余弦距离Cosine distance(‘cosine’)与Jaccard距离相比，Cosine距离不仅忽略0-0匹配，而且能够处理非二元向量，即考虑到变量值的大小。
    %对这两者，距离与相似度和为一。
    %https://www.cnblogs.com/chaosimple/archive/2013/06/28/3160839.html

    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
%     Distance = pdist2(PopObj,Z,'euclidean');%上面注释的两个语句也可以哦
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);%除了最后front的每一个参考点，计算关联解的个数

    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));%所有参考点中与之关联个数最少的那个参考点
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end