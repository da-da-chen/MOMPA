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
function [fobj, fcon,P] = mompa_getMOFcn(F, decisionVar, numObj)

[N,D]  = size(decisionVar);
switch F
    case 'ZDT1'
        fobj(:, 1) = decisionVar(:, 1);
        g = 1 + 9*mean(decisionVar(:,2:end),2);
        h = 1 - (fobj(:,1)./g).^0.5;
        fobj(:,2) = g.*h;

        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;
    case 'ZDT2'
        fobj(:,1) = decisionVar(:,1);
        g = 1 + 9*mean(decisionVar(:,2:end),2);
        h = 1 - (fobj(:,1)./g).^2;
        fobj(:,2) = g.*h;
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^2;
    case 'ZDT3'
        fobj(:,1) = decisionVar(:,1);
        g = 1 + 9*mean(decisionVar(:,2:end),2);
        h = 1 - (fobj(:,1)./g).^0.5 - fobj(:,1)./g.*sin(10*pi*fobj(:,1));
        fobj(:,2) = g.*h;
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5 - P(:,1).*sin(10*pi*P(:,1));
        P      = P(NDSort(P,fcon,1)==1,:);
    case 'ZDT4'
        fobj(:,1) = decisionVar(:,1);
        g = 1 + 10*(size(decisionVar,2)-1) + sum(decisionVar(:,2:end).^2-10*cos(4*pi*decisionVar(:,2:end)),2);
        h = 1 - (fobj(:,1)./g).^0.5;
        fobj(:,2) = g.*h;
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P(:,1) = (0:1/(N-1):1)';
        P(:,2) = 1 - P(:,1).^0.5;

    case 'ZDT6'
        fobj(:,1) = 1 - exp(-4*decisionVar(:,1)).*sin(6*pi*decisionVar(:,1)).^6;
        g = 1 + 9*mean(decisionVar(:,2:end),2).^0.25;
        h = 1 - (fobj(:,1)./g).^2;
        fobj(:,2) = g.*h;
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        minf1  = 0.280775;
        P(:,1) = (minf1:(1-minf1)/(N-1):1)';
        P(:,2) = 1 - P(:,1).^2;

    case 'DTLZ1'
        M      = numObj;
        g      = 100*(D-M+1+sum((decisionVar(:,M:end)-0.5).^2-cos(20.*pi.*(decisionVar(:,M:end)-0.5)),2));
        fobj = 0.5*repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),decisionVar(:,1:M-1)],2)).*[ones(N,1),1-decisionVar(:,M-1:-1:1)];
        %......................
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P = UniformPoint(N,numObj)/2;
    case 'DTLZ2'
        M = numObj;
        g = sum((decisionVar(:,M:end)-0.5).^2,2);
        fobj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(decisionVar(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(decisionVar(:,M-1:-1:1)*pi/2)];
        %......................
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P = UniformPoint(N,numObj);
        P = P./repmat(sqrt(sum(P.^2,2)),1,numObj);
    case 'DTLZ3'
        M      = numObj;
        g      = 100*(D-M+1+sum((decisionVar(:,M:end)-0.5).^2-cos(20.*pi.*(decisionVar(:,M:end)-0.5)),2));
        fobj = repmat(1+g,1,M).*fliplr(cumprod([ones(N,1),cos(decisionVar(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(decisionVar(:,M-1:-1:1)*pi/2)];
        %......................
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
    case 'DTLZ4'
        M      = numObj;
        decisionVar(:,1:M-1) = decisionVar(:,1:M-1).^100;
        g      = sum((decisionVar(:,M:end)-0.5).^2,2);
        fobj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(decisionVar(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(decisionVar(:,M-1:-1:1)*pi/2)];
        %......................
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
    case 'DTLZ5'
        M      = numObj;
        g      = sum((decisionVar(:,M:end)-0.5).^2,2);
        Temp   = repmat(g,1,M-2);
        decisionVar(:,2:M-1) = (1+2*Temp.*decisionVar(:,2:M-1))./(2+2*Temp);
        fobj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(decisionVar(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(decisionVar(:,M-1:-1:1)*pi/2)];
        %......................
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P = [0:1/(N-1):1;1:-1/(N-1):0]';
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-2)),P];
        P = P./sqrt(2).^repmat([M-2,M-2:-1:0],size(P,1),1);
    case 'DTLZ6'
        M      = numObj;
        g      = sum(decisionVar(:,M:end).^0.1,2);
        Temp   = repmat(g,1,M-2);
        decisionVar(:,2:M-1) = (1+2*Temp.*decisionVar(:,2:M-1))./(2+2*Temp);
        fobj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(decisionVar(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(decisionVar(:,M-1:-1:1)*pi/2)];
        %......................
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        P = [0:1/(N-1):1;1:-1/(N-1):0]';
        P = P./repmat(sqrt(sum(P.^2,2)),1,size(P,2));
        P = [P(:,ones(1,M-2)),P];
        P = P./sqrt(2).^repmat([M-2,M-2:-1:0],size(P,1),1);
    case 'DTLZ7'
        M               = numObj;
        fobj         = zeros(size(decisionVar,1),M);
        g               = 1+9*mean(decisionVar(:,M:end),2);
        fobj(:,1:M-1) = decisionVar(:,1:M-1);
        fobj(:,M)     = (1+g).*(M-sum(fobj(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*fobj(:,1:M-1))),2));
        %......................
        fcon(:, 1:D) = -decisionVar;
        fcon(:, D+1: 2*D) = decisionVar - 1;

        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        X            = ReplicatePoint(N,M-1);
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        P            = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
    case 'WFG1'
        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = WFG_s_linear(z01(:,K+1:end),0.35);

        t2 = zeros(N,K+L);
        t2(:,1:K)     = t1(:,1:K);
        t2(:,K+1:end) = WFG1_b_flat(t1(:,K+1:end),0.8,0.75,0.85);

        t3 = zeros(N,K+L);
        t3 = WFG1_b_poly(t2,0.02);

        t4 = zeros(N,M);
        for i = 1 : M-1
            t4(:,i) = WFG_r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)),2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
        end
        t4(:,M) = WFG_r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t4(:,M),A(i)).*(t4(:,i)-0.5)+0.5;
        end
        x(:,M) = t4(:,M);

        h      = WFG_convex(x);
        h(:,M) = WFG1_mixed(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        %......................
        fcon = -decisionVar;


        P = UniformPoint(N,M);
        c = ones(size(P,1),M);
        for i = 1 : size(P,1)
            for j = 2 : M
                temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a+cos(10*pi*a+pi/2)/10/pi,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        P      = WFG_convex(x);
        P(:,M) = WFG1_mixed(x);
        P      = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG2'
        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = WFG_s_linear(z01(:,K+1:end),0.35);

        t2 = zeros(N,K+L/2);
        t2(:,1:K) = t1(:,1:K);
        t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;

        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = WFG_r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = WFG_r_sum(t2(:,K+1:K+L/2),ones(1,L/2));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);

        h      = WFG_convex(x);
        h(:,M) = WFG2_disc(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        %......................
        fcon = -decisionVar;

        P = UniformPoint(N,M);
        c = ones(size(P,1),M);
        for i = 1 : size(P,1)
            for j = 2 : M
                temp = P(i,j)/P(i,1)*prod(1-c(i,M-j+2:M-1));
                c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*P(:,M)./P(:,M-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a.*cos(5*pi*a).^2,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        P      = WFG_convex(x);
        P(:,M) = WFG2_disc(x);
        P      = P(NDSort(P,fcon,1)==1,:);
        P      = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG3'
        %......................
        fcon = -decisionVar;

        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = [1,zeros(1,M-2)];

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = WFG_s_linear(z01(:,K+1:end),0.35);

        t2 = zeros(N,K+L/2);
        t2(:,1:K) = t1(:,1:K);
        t2(:,K+1:K+L/2) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;

        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = WFG_r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = WFG_r_sum(t2(:,K+1:K+L/2),ones(1,L/2));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);

        h      = WFG3_linear(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        X = (0:1/(N-1):1)';
        X = [X,zeros(N,M-2)+0.5,zeros(N,1)];
        P = WFG3_linear(X);
        P = repmat(2:2:2*M,size(P,1),1).*P;

    case 'WFG4'
        %......................
        fcon = -decisionVar;

        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        t1 = WFG4_s_multi(z01,30,10,0.35);

        t2 = zeros(N,M);
        for i = 1 : M-1
            t2(:,i) = WFG_r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t2(:,M) = WFG_r_sum(t1(:,K+1:K+L),ones(1,L));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
        end
        x(:,M) = t2(:,M);

        h = WFG4_concave(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG5'
        %......................
        fcon = -decisionVar;

        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        t1 = WFG5_s_decept(z01,0.35,0.001,0.05);

        t2 = zeros(N,M);
        for i = 1 : M-1
            t2(:,i) = WFG5_r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t2(:,M) = WFG5_r_sum(t1(:,K+1:K+L),ones(1,L));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
        end
        x(:,M) = t2(:,M);

        h = WFG5_concave(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG6'
        %......................
        fcon = -decisionVar;

        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        t1(:,1:K)     = z01(:,1:K);
        t1(:,K+1:end) = WFG6_s_linear(z01(:,K+1:end),0.35);

        t2 = zeros(N,M);
        for i = 1 : M-1
            t2(:,i) = WFG6_r_nonsep(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
        end

        SUM = zeros(N,1);
        for i = K+1 : K+L-1
            for j = i+1 : K+L
                SUM = SUM + abs(t1(:,i)-t1(:,j));
            end
        end
        t2(:,M) = (sum(t1(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
        end
        x(:,M) = t2(:,M);

        h = WFG6_concave(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG7'
        %......................
        fcon = -decisionVar;

        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
        t1(:,1:K) = z01(:,1:K).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K)).*abs(floor(0.5-Y(:,1:K))+0.98/49.98)));
        t1(:,K+1:end) = z01(:,K+1:end);

        t2 = zeros(N,K+L);
        t2(:,1:K)     = t1(:,1:K);
        t2(:,K+1:end) = WFG7_s_linear(t1(:,K+1:end),0.35);

        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = WFG7_r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = WFG7_r_sum(t2(:,K+1:K+L),ones(1,L));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);

        h = WFG7_concave(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG8'
        %......................
        fcon = -decisionVar;

        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        t1(:,1:K) = z01(:,1:K);
        Y = (cumsum(z01,2)-z01)./repmat(0:K+L-1,N,1);
        t1(:,K+1:K+L) = z01(:,K+1:K+L).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,K+1:K+L)).*abs(floor(0.5-Y(:,K+1:K+L))+0.98/49.98)));

        t2 = zeros(N,K+L);
        t2(:,1:K)     = t1(:,1:K);
        t2(:,K+1:end) = WFG8_s_linear(t1(:,K+1:end),0.35);

        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = WFG8_r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
        end
        t3(:,M) = WFG8_r_sum(t2(:,K+1:K+L),ones(1,L));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);

        h = WFG8_concave(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = repmat(2:2:2*M,size(P,1),1).*P;
    case 'WFG9'
        %......................
        fcon = -decisionVar;

        M               = numObj;
        K = M-1;
        L = D - K;
        D = 1;
        S = 2 : 2 : 2*M;
        A = ones(1,M-1);

        z01 = decisionVar./repmat(2:2:size(decisionVar,2)*2,N,1);

        t1 = zeros(N,K+L);
        Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
        t1(:,1:K+L-1) = z01(:,1:K+L-1).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K+L-1)).*abs(floor(0.5-Y(:,1:K+L-1))+0.98/49.98)));
        t1(:,end)     = z01(:,end);

        t2 = zeros(N,K+L);
        t2(:,1:K)     = WFG9_s_decept(t1(:,1:K),0.35,0.001,0.05);
        t2(:,K+1:end) = WFG9_s_multi(t1(:,K+1:end),30,95,0.35);

        t3 = zeros(N,M);
        for i = 1 : M-1
            t3(:,i) = WFG9_r_nonsep(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
        end
        SUM = zeros(N,1);
        for i = K+1 : K+L-1
            for j = i+1 : K+L
                SUM = SUM + abs(t2(:,i)-t2(:,j));
            end
        end
        t3(:,M) = (sum(t2(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));

        x = zeros(N,M);
        for i = 1 : M-1
            x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
        end
        x(:,M) = t3(:,M);

        h = WFG9_concave(x);
        fobj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

        P = UniformPoint(N,M);
        P = P./repmat(sqrt(sum(P.^2,2)),1,M);
        P = repmat(2:2:2*M,size(P,1),1).*P;

end
end

function W = ReplicatePoint(SampleNum,M)
if M > 1
    SampleNum = (ceil(SampleNum^(1/M)))^M;
    Gap       = 0:1/(SampleNum^(1/M)-1):1;
    eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
    eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
else
    W = (0:1/(SampleNum-1):1)';
end
end
function Output = WFG_s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = WFG1_b_flat(y,A,B,C)
Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
Output = roundn(Output,-6);
end

function Output = WFG1_b_poly(y,a)
Output = y.^a;
end

function Output = WFG_r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = WFG_convex(x)
Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = WFG1_mixed(x)
Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end
function Output = WFG2_disc(x)
Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end
function Output = WFG3_linear(x)
Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end
function Output = WFG4_s_multi(y,A,B,C)
Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end
function Output = WFG4_concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
function Output = WFG5_s_decept(y,A,B,C)
Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end

function Output = WFG5_r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = WFG5_concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
function Output = WFG6_s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = WFG6_r_nonsep(y,A)
Output = zeros(size(y,1),1);
for j = 1 : size(y,2)
    Temp = zeros(size(y,1),1);
    for k = 0 : A-2
        Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
    end
    Output = Output+y(:,j)+Temp;
end
Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end

function Output = WFG6_concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
function Output = WFG7_s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = WFG7_r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = WFG7_concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
function Output = WFG8_s_linear(y,A)
Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = WFG8_r_sum(y,w)
Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = WFG8_concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
function Output = WFG9_s_decept(y,A,B,C)
Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end

function Output = WFG9_s_multi(y,A,B,C)
Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end

function Output = WFG9_r_nonsep(y,A)
Output = zeros(size(y,1),1);
for j = 1 : size(y,2)
    Temp = zeros(size(y,1),1);
    for k = 0 : A-2
        Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
    end
    Output = Output+y(:,j)+Temp;
end
Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end

function Output = WFG9_concave(x)
Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
function W = VNT_ReplicatePoint(SampleNum,M)
if M > 1
    SampleNum = (ceil(SampleNum^(1/M)))^M;
    Gap       = 0:1/(SampleNum^(1/M)-1):1;
    eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
    eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
else
    W = (0:1/(SampleNum-1):1)';
end
end