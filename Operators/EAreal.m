function OffspringDec = EAreal(ParentDec,D,lu)

%%%Parameters
proC=0.9;
disC=10;
proM=1;
disM=20;

[N,ND] = size(ParentDec);

Parent1Dec = ParentDec(1:N/2,1:D);
Parent2Dec = ParentDec(N/2+1:end,1:D);
 %% Simulated binary crossover
beta = zeros(N/2,D);
mu   = rand(N/2,D);
beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
beta = beta.*(-1).^randi([0,1],N/2,D);
beta(rand(N/2,D)<0.5) = 1;
beta(repmat(rand(N/2,1)>proC,1,D)) = 1;
OffspringDec = [(Parent1Dec+Parent2Dec)/2+beta.*(Parent1Dec-Parent2Dec)/2
    (Parent1Dec+Parent2Dec)/2-beta.*(Parent1Dec-Parent2Dec)/2];


 %% Polynomial mutation
    Lower = repmat(lu(1,:),N,1);
    Upper = repmat(lu(2,:),N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                         (1-(OffspringDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    OffspringDec(temp) = OffspringDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                         (1-(Upper(temp)-OffspringDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
                                     