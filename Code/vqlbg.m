function [m, p, DistHist]=vqlbg(X,L)

% X: a matrix each column of which is a data vector

e=.01;     % X---> [X-e*X and X+e*X] Percentage for Spliting
eRed=0.75; % e=e*eRed;
DT=.005; % The threshold in improvement in Distortion before terminating and spliting again
DTRed=0.75; % Rate of reduction of Improvement Threshold, DT, after each spliting
MinPop=0.10; 
d=size(X,1); % Dimension
N=size(X,2); % Number of Data points
isFirstRound=1; % First Iteration after Spliting
if numel(L)==1
    M=mean(X,2); % Mean Vector
    CB=[M*(1+e) M*(1-e)]; % Split to two vectors
else
    CB=L; % If the codebook is passed to the function just train it
    L=size(CB,2);
    e=e*(eRed^fix(log2(L)));
    DT=DT*(DTRed^fix(log2(L)));
end
LC=size(CB,2); % Current size of the codebook
Iter=0;
Split=0;
IsThereABestCB=0;
maxIterInEachSize=20;
EachSizeIterCounter=0;
while 1
    %Distance Calculation
    [minIndx, dst]=VQIndex(X,CB); % Find the closest codewords to each data vector
    ClusterD=zeros(1,LC);
    Population=zeros(1,LC);
    LowPop=[];
    % Find the Centroids (Mean of each Cluster)
    for i=1:LC
        Ind=find(minIndx==i);
        if length(Ind)<MinPop*N/LC %if a cluster has very low population just remember it
            LowPop=[LowPop i];
        else
            CB(:,i)=mean(X(:,Ind),2);
            Population(i)=length(Ind);
            ClusterD(i)=sum(dst(Ind));
        end        
    end
    if ~isempty(LowPop)
        [temp MaxInd]=maxn(Population,length(LowPop));
        CB(:,LowPop)=CB(:,MaxInd)*(1+e); 
        CB(:,MaxInd)=CB(:,MaxInd)*(1-e);
        
        [minIndx, dst]=VQIndex(X,CB);
        ClusterD=zeros(1,LC);
        Population=zeros(1,LC);
        
        for i=1:LC
            Ind=find(minIndx==i);
            if ~isempty(Ind)
                CB(:,i)=mean(X(:,Ind),2);
                Population(i)=length(Ind);
                ClusterD(i)=sum(dst(Ind));
            else 
                CB(:,i)=X(:,fix(rand*N)+1);
                disp('A random vector was assigned as a codeword.')
                isFirstRound=1;% At least another iteration is required
            end                
        end
    end
    Iter=Iter+1;
    if isFirstRound % First iteration after a split (dont exit)
        TotalDist=sum(ClusterD(~isnan(ClusterD)));
        DistHist(Iter)=TotalDist;
        PrevTotalDist=TotalDist;        
        isFirstRound=0;
    else
        TotalDist=sum(ClusterD(~isnan(ClusterD)));  
        DistHist(Iter)=TotalDist;
        PercentageImprovement=((PrevTotalDist-TotalDist)/PrevTotalDist);
        if PercentageImprovement>=DT 
            PrevTotalDist=TotalDist; 
            isFirstRound=0;
        else
            EachSizeIterCounter=0;
            if LC>=L 
                if L==LC %Exact number of codewords
                    disp(TotalDist)
                    break
                else % Remove one codeword at a time
                    [temp, Ind]=min(Population); % Eliminate low population codewords
                    NCB=zeros(d,LC-1);
                    NCB=CB(:,setxor(1:LC,Ind(1)));
                    CB=NCB;
                    LC=LC-1;
                    isFirstRound=1;
                end
            else %If not enough codewords yet, then Split more
                CB=[CB*(1+e) CB*(1-e)];
                e=eRed*e; %Split size reduction
                DT=DT*DTRed; %Improvement Threshold Reduction
                LC=size(CB,2);
                isFirstRound=1;
                Split=Split+1;
                IsThereABestCB=0;
                disp(LC)
            end
        end
    end    
    if ~IsThereABestCB
        BestCB=CB;
        BestD=TotalDist;
        IsThereABestCB=1;
    else
        if TotalDist<BestD
            BestCB=CB;
            BestD=TotalDist;
        end
    end
    EachSizeIterCounter=EachSizeIterCounter+1;
    if EachSizeIterCounter>maxIterInEachSize 
        EachSizeIterCounter=0;
        CB=BestCB; 
        IsThereABestCB=0;
        if LC>=L 
            if L==LC 
                disp(TotalDist)
                break
            else 
                [temp, Ind]=min(Population);
                NCB=zeros(d,LC-1);
                NCB=CB(:,setxor(1:LC,Ind(1)));
                CB=NCB;
                LC=LC-1;
                isFirstRound=1;
            end
        else
            CB=[CB*(1+e) CB*(1-e)];
            e=eRed*e;
            DT=DT*DTRed; 
            LC=size(CB,2);
            isFirstRound=1;
            Split=Split+1;
            IsThereABestCB=0;
            disp(LC)
        end
    end        
    disp(TotalDist)
    p=Population/N;
    save CBTemp CB p DistHist
end
m=CB;
p=Population/N;
disp(['Iterations = ' num2str(Iter)])
disp(['Split = ' num2str(Split)])
function [v, i]=maxn(x,n)
if nargin<2
    [v, i]=max(x);
else
    n=min(length(x),n);
    [v, i]=sort(x);
    v=v(end:-1:end-n+1);
    i=i(end:-1:end-n+1);    
end
        
function [I, dst]=VQIndex(X,CB) 
L=size(CB,2);
N=size(X,2);
LNThreshold=64*10000;
if L*N<LNThreshold
    D=zeros(L,N);
    for i=1:L
        D(i,:)=sum((repmat(CB(:,i),1,N)-X).^2,1);
    end
    [dst I]=min(D);
else
    I=zeros(1,N);
    dst=I;
    for i=1:N
        D=sum((repmat(X(:,i),1,L)-CB).^2,1);
        [dst(i) I(i)]=min(D);
    end
end
    
function [I, dist]=VQLSFSpectralIndex(X,CB,W)
if nargin<3
    L=256;
    W=ones(L,1);
else
    if isscalar(W)
        L=W;
        W=ones(L,1);
    elseif isvector(W)
        W=W(:);
        L=length(W);
    else
        error('Invalid input argument. W should be either a vector or a scaler!')
    end
end
NX=size(X,2);
NCB=size(CB,2);
AX=lsf2lpc(X);
ACB=lsf2lpc(CB);
D=zeros(NCB,1);
w=linspace(0,pi,L+1);
w=w(1:end-1);
N=size(AX,2)-1;
WFZ=zeros(N+1,L);
IMAGUNIT=sqrt(-1);
for k=0:N
    WFZ(k+1,:)=exp(IMAGUNIT*k*w);
end
SCB=zeros(L,NCB);
for i=1:NCB
    SCB(:,i)=(1./abs(ACB(i,:)*WFZ));
end
I=zeros(1,NX);
dist=zeros(1,NX);
for j=1:NX
    SX=(1./abs(AX(j,:)*WFZ))';    
    for i=1:NCB
        D(i)=sqrt(sum(((SX-SCB(:,i)).^2).*W));
    end
    [dist(j), I(j)]=min(D);
end