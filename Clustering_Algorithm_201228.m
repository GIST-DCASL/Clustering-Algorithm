clear all;
clc;
%% Read file
% Input is the incidence matrix of the graph.
% Prepare the Excel file of incidence matrix which all coefficents are 1.
% (1,1) component of the incidence matrix is placed at A1.

filename = 'Example1';
Inc = xlsread(filename,'A1:ZZ30');
[node,edge] = size(Inc);

%% 1. No '1' in the row of incidence matrix -> leader node
% Class-1 node is the node which has no in-degree edge.
% These kinds of nodes are leader node and it is a CR node.
% Check whether each row of incidence matrix does not have 1-component.
% Step 1 of Clustering Algorithm

count1 = 0;
R1N = zeros(1,node);
for i = 1 : node
    Row = Inc(i,:);
    CheckRow = ( Row == 1 );
    if ( CheckRow == zeros(1,edge) )        % Classify class-1 nodes
       count1 = count1 + 1;
       R1N(1,count1) = i;
    end
end
r1 = find(R1N);
R1node = R1N(1:numel(r1));

%% 2.Only one '1' in the incidence matrix -> follower node
% Class-2 node is the node which has only one in-degree edge.
% These kinds of nodes are follower node and it is a CF node.
% Check whether each row of incidence matrix has only one 1-component.
% Step 2 of Clustering Algorithm

count2 = 0;
R2N = zeros(1,node);
for j = 1 : node
    Row = Inc(j,:);
    CheckRow = ( Row == 1 );
    S = sum(CheckRow);
    if (S == 1)                             % Classify class-2 nodes
        count2 = count2 + 1;
        R2N(1,count2) = j;
    end
end
r2 = find(R2N);
R2node = R2N(1:numel(r2));

%% 3. Several(>=2) '1' in the incidence matrix
% Class-3 node is the node which has more than or equal to two in-degree edges.
% These kinds of nodes can become CR or CF node.
% So I desgin and use R3fnc function
% Step 3 & 4 & 5 of Clustering Algorithm


count3 = 0;
R3N = zeros(1,node);
for k = 1 : node
    Row = Inc(k,:);
    CheckRow = ( Row == 1 );
    S = sum(CheckRow);
    if (S >= 2)                             %Classifiy class-3 nodes
        count3 = count3 + 1;
        R3N(1,count3) = k;
    end
end

R3R = zeros(1,node);
R3F = zeros(1,node);
numfollower = 0;
numroot = 0;

for t = 1:node
    if( sum(R3N == t) ~= 0 )
        Wlast = [t];
        [X,Y] = R3fnc(Inc,t,R3N,edge,node,R1node,R2node,Wlast);         % Traverse all possible non-repeated reverse directed path
        Len = length(X);
        [R,C] = size(Y);
        if ( sum( X == X(1,1) ) == Len )
            numfollower = numfollower + 1;
            R3F(1,numfollower) = t;
        else
            check = 0;
            for i = 1:node
                if ( i ~= t )
                    for j = 1:R
                        if ( sum(Y(j,:) == i) >= 1 )        %Check common passed node
                            check = check + 1;
                        end
                    end
                    if (check ~= R)
                        check = 0;
                    else
                        break;
                    end
                end
            end
            if ( check == R )                           % Collect CF node in class-3
                numfollower = numfollower + 1;
                R3F(1,numfollower) = t;
            else                                        % Collect CR node in class-3
                numroot = numroot + 1;
                R3R(1,numroot) = t;
            end
        end
    end
end

if ( R3R == zeros(1,node) )
    R3root = [0];
else
    r3r = find(R3R);
    R3root = R3R(1:numel(r3r));
end

if ( R3F == zeros(1,node) )
    R3follow = [0];
else
    r3f = find(R3F);
    R3follow = R3F(1:numel(r3f));
end



%% Clustereing
% Print out the clustering result.
% Step 6 of Clustering Algorithm

Clusternum = 0;
Clusterset = zeros(node, node);

for nodenum = 1 : node
    if ( (sum(R1node == nodenum) + sum(R3root == nodenum) ) ~= 0 )
        Clusternum = Clusternum + 1;                                    % Count number of cluster
        Clusterset(nodenum,1) = nodenum;                                % Take a CR node at each cluster
    end
end

for nodenum = 1 : node
    if ( sum(R2node == nodenum) ~= 0 )
        numcnt = 1;
        for e = 1 : edge
            if ( Inc(nodenum,e) == 1 )
                for n = 1 : node
                    if ( Inc(n,e) == -1 )
                        if ( ( sum(R2node == n) + sum(R3follow == n) ) ~= 0 )
                            Wlast = [n];
                            X = follower(Inc,n,edge,node,R1node,R2node,R3follow,R3root,Wlast);
                            n = X(1,1);
                        end
                        for m = 1 : node
                            if ( Clusterset(n,m) ~= 0 )
                                numcnt = numcnt + 1 ;
                            end
                        end    
                        Clusterset(n,numcnt) = nodenum ;                % Match a class-2 CF node with appropriate CR node.
                    end
                end
            end
        end
    end
end

if( sum(R3follow) ~= 0 )
    Len = length(R3follow);
    for r = 1 : Len
        k = R3follow(1,r);
        Wlast = [k];
        RNode = follower(Inc,k,edge,node,R1node,R2node,R3follow,R3root,Wlast);
        RNod = RNode(1,1);
        numcnt = 1;
        for m = 1 : node
            if ( Clusterset(RNod,m) ~= 0 )
                numcnt = numcnt + 1 ;
            end
        end
        Clusterset(RNod,numcnt) = k ;                                   % Match a class-3 CF node with appropriate CR node.
    end
end


Clusternum
Clusterset


%% R3node function
% This function is a recursive function to traverse all possible non-repeated reverse directed path.
% This function gives first-reached CR nodes & all nodes on the reverse directed path.

function [X,Y] = R3fnc(Inc,t,R3N,edge,node,R1node,R2node,Wlast)
if( sum(R3N == t) ~= 0 )
    kcnt = 1;
    pcnt = 1;
    X=0;
    Y=0;
    for e = 1 : edge
        if ( Inc(t,e) == 1 )
            for n = 1 : node
                if ( Inc(n,e) == -1 )
                   if ( sum(Wlast == n) == 0 )
                        if ( sum(R1node == n) ~= 0 )
                            X(1,kcnt) = n;
                            Y(pcnt,1) = t;
                            Y(pcnt,2) = n;
                            kcnt = kcnt + 1;
                            pcnt = pcnt + 1;
                        else
                            N(1,1) = n;
                            Wlast = horzcat(Wlast,N);
                            [L, W] = R3fnc(Inc,n,R3N,edge,node,R1node,R2node,Wlast);
                            if ( sum(L,'all')~=0 && sum(W,'all')~=0 ) 
                                w = length(L);
                                for b = 1:w
                                    X(1,kcnt+b-1) = L(1,b);
                                end
                                kcnt = kcnt + w;
                                [p,q] = size(W);
                                for c = 1:p
                                    for d = 1:q
                                        Y(c+pcnt-1,1+d) = W(c,d);
                                    end
                                    Y(c+pcnt-1,1) = t;
                                end
                                pcnt = pcnt + p;
                            end
                            LenW = length(Wlast);
                            Wlast(:,LenW)=[];
                        end
                   else

                   end
                end
            end
        end
    end
end
if ( sum(R2node == t) ~= 0 )
    kcnt = 1;
    pcnt = 1;
    X=0;
    Y=0;
    for e = 1 : edge
        if ( Inc(t,e) == 1 )
            for n = 1 : node
                if ( Inc(n,e) == -1 )
                    if ( sum(Wlast == n) == 0 )
                        if ( sum(R1node == n) ~= 0 )
                            X(1,kcnt) = n;
                            Y(kcnt,1) = t;
                            Y(kcnt,2) = n;
                            kcnt = kcnt + 1;
                            pcnt = pcnt + 1;
                        else
                            N(1,1) = n;
                            Wlast = horzcat(Wlast,N);
                            [L, W] = R3fnc(Inc,n,R3N,edge,node,R1node,R2node,Wlast);
                            w = length(L);
                            for b = 1:w
                                X(1,kcnt+b-1) = L(1,b);
                            end
                            kcnt = kcnt + w;
                            [p,q] = size(W);
                            for c = 1:p
                                for d = 1:q
                                    Y(c+pcnt-1,1+d) = W(c,d);
                                end
                                Y(c+pcnt-1,1) = t;
                            end
                            pcnt = pcnt + p;
                            LenW = length(Wlast);
                            Wlast(:,LenW)=[];
                        end
                    end
                    
                end
            end
        end
    end
end
end


%% Clustering function
% This funciton helps us to match the class-3 CF node with appropriate CR node.
% This function is a recursive function to find appropriate CR node for each CF node.

function X = follower(Inc,t,edge,node,R1node,R2node,R3follow,R3root,Wlast)
if( ( sum(R2node == t) + sum(R3follow == t) ) ~= 0 )
    kcnt = 1;
    X=0;
    for e = 1 : edge
        if ( Inc(t,e) == 1 )
            for n = 1 : node
                if ( Inc(n,e) == -1 )
                   if ( sum(Wlast == n) == 0 )
                        if ( (sum(R1node == n) + sum(R3root == n) ) ~= 0 )
                            X(1,kcnt) = n;
                            kcnt = kcnt + 1;
                        else
                            N(1,1) = n;
                            Wlast = horzcat(Wlast,N);
                            L = follower(Inc,n,edge,node,R1node,R2node,R3follow,R3root,Wlast);
                            if ( sum(L,'all')~=0 )
                                w = length(L);
                                for b = 1:w
                                    X(1,kcnt+b-1) = L(1,b);
                                end
                                kcnt = kcnt + w;
                            end
                            LenW = length(Wlast);
                            Wlast(:,LenW)=[];
                        end
                   else
                       
                   end
                end
            end
        end
    end
end
end