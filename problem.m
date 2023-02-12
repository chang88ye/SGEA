function prob=problem(prob_id)
%% Implementation of popular test suites
% test suites are numbered, i.e. CEC2018-DF (1), FDA (2), dCOEA-DMOP (3), PPS-F (4), etc
%
% CEC2018-DF: ID 101 - 114
%             - The 14 test functions DF1 - DF14 are for cec2018 competition on dynamic multiobjective optimisation. 
% 
% FDA      :  ID 201 - 205
%             - The 5 test functions FDA1 - FDA5 pioneered research on EDMO
% 
% DMOP     :  ID 301 - 303
%             - three test functions DMOP1 - DMOP3 in dCOEA
%           
% F        :  ID 401 - 408
%             - Eight test functions F1 - F8 in PPS
    
    prob=struct('dynPar', [100,  10,  10]); % create a struct and provide values for TO, taut, and nt

    switch prob_id

        %% DF test suite: CEC2018
        case 101                 % 'DF1' (dMOP2)
            prob.func=@DF1;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 102                 % DF2 (modified dMOP3)
            prob.func=@DF2;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 103                 % 'DF3' (ZJZ)
            
            prob.func=@DF3;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 2*ones(1,prob.varDim-1)];
        case 104                 % DF4
            
            prob.func=@DF4;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[-2*ones(1,prob.varDim);2*ones(1,prob.varDim)];
        case 105                 % DF5 (modified JY2)
            
            prob.func=@DF5;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 106                 % DF6  (modified JY7)
            
            prob.func=@DF6;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 107                 % DF7
            
            prob.func=@DF7;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[1 zeros(1,prob.varDim-1);4 ones(1,prob.varDim-1)];
        case 108                 % DF8
            
            prob.func=@DF8;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 109            % DF9
            
            prob.func=@DF9;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 110   %DF10
            
            prob.func=@DF10;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        case 111   %DF11
            
            prob.func=@DF11;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 112 % DF12
            
            prob.func=@DF12;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        case 113 %DF13 ->time consuming
           
            prob.func=@DF13;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        case 114 % DF14
           
            prob.func=@DF14;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
    
        %% FDA test suite: 
        case 201 % FDA1
            
            prob.func=@FDA1;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 202 %  FDA2   将xiii看成空集
            
            prob.func=@FDA2;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[0 -ones(1,prob.varDim-1);1 ones(1,prob.varDim-1)];
        case 203 %  FDA3
            
            prob.func=@FDA3;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) ones(1,prob.varDim-2)];
        
        case 204 % FDA4
            
            prob.func=@FDA4;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 205   % FDA5

            prob.func=@FDA5;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];

        %% DMOP test suite: dCOEA
        case 301 % dmop1
            prob.func=@DMOP1;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];

        case 302 %%%%dmop2
            
            prob.func=@DMOP2;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
        case 303    % dMOP3
            
            prob.func=@DMOP3;
            prob.varDim=10;
            prob.objDim=2;
            prob.bounds=[zeros(1,prob.varDim);ones(1,prob.varDim)];
            
       %% F test suite: A population prediction strategy for evolutionary dynamic multiobjective optimization
       case 408   %%%F8
            
            prob.func=@F8;
            prob.varDim=10;
            prob.objDim=3;
            prob.bounds=[zeros(1,2) -ones(1,prob.varDim-2);ones(1,2) 2*ones(1,prob.varDim-2)];

    end
    prob.update_time=@cal_time;
end




%% function to calculate time instant
function t=cal_time(tau, taut, nt, T0)
% --------------------------------------------------------|%
% The "time" term in the test suite is defined as:        |
%          t=1/nt*floor(tau/taut)                         |
% where - nt:    severity of change                       |
%       - taut:  frequency of change                      |
%       - tau:   current generation counter               |
% --------------------------------------------------------|%
% the first change occurs after T0 generations, that is, the
%generation at which a change occurs is (T0+1), (T0+taut+1), etc.  

    tau_tmp=max(tau+taut-(T0+1),0);
    t=1/nt*floor(tau_tmp/taut);
end


%% test problem definition
function xf= DF1(x,t)
    f=[];
    G=abs(sin(0.5*pi*t));
    H=0.75*sin(0.5*pi*t)+1.25;
    g=1+sum((x(:,2:end)-G).^2,2);
    f1=x(:,1);
    f2=g.*(1-(x(:,1)./g).^H);
    f=[f1,f2];
    xf=[x,f];
end

function xf=DF2(x,t)
    f=[];
    G=abs(sin(0.5*pi*t));
    r=1+floor((n-1)*G);
    tmp=setdiff(1:n,r);
    g=1+sum((x(:,tmp)-G).^2,2);
    f1=x(:,r);
    f2=g.*(1-(x(:,r)./g).^0.5);
    f=[f1,f2];
end

function f= DF3(x,t)
    f=[];
    G=sin(0.5*pi*t);
    H=G+1.5;
    g=1+sum((x(:,2:end)-G-x(:,1).^H).^2,2);
    f1=x(:,1);
    f2=g.*(1-(x(:,1)./g).^H);
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF4(x,t)
    f=[];
    n=size(x,2);

    a=sin(0.5*pi*t);
    b=1+abs(cos(0.5*pi*t));
    c=max(abs(a), a+b);
    H=1.5+a;
    g=1+sum((x(:,2:end)-a*(x(:,1)/c).^2./[2:n]).^2,2);
    f1=g.*abs(x(:,1)-a).^H;
    f2=g.*abs(x(:,1)-a-b).^H;
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF5(x,t)
    f=[];
    G=sin(0.5*pi*t);
    w=floor(10*G);
    g=1+sum((x(:,2:end)-G).^2,2);
    f1=g.*(x(:,1)+0.02*sin(w*pi*x(:,1)));
    f2=g.*(1-x(:,1)+0.02*sin(w*pi*x(:,1)));
    f=[f1,f2];
    xf=[x,f];
end

function f= DF6(x,t)
    f=[];
    G=sin(0.5*pi*t);
    a=0.2+2.8*abs(G);
    y=x(:,2:end)-G;
    g=1+sum((abs(G)*y.^2-10*cos(2*pi*y)+10),2);
    f1=g.*(x(:,1)+0.1*sin(3*pi*x(:,1))).^a;
    f2=g.*(1-x(:,1)+0.1*sin(3*pi*x(:,1))).^a;
    f=[f1,f2];
end

function xf= DF7(x,t)
    f=[];
    a=5*cos(0.5*pi*t);
    tmp=1./(1+exp(a*(x(:,1)-2.5)));
    g=1+sum((x(:,2:end)-tmp).^2,2);
    f1=g.*(1+t)./x(:,1);
    f2=g.*x(:,1)/(1+t);
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF8(x,t)
    f=[];
    G=sin(0.5*pi*t);
    a=2.25+2*cos(2*pi*t);
    b=100*G^2;
    tmp=G*sin(4*pi*x(:,1).^b)/(1+abs(G));
    g=1+sum((x(:,2:end)-tmp).^2,2);
    f1=g.*(x(:,1)+0.1*sin(3*pi*x(:,1)));
    f2=g.*(1-x(:,1)+0.1*sin(3*pi*x(:,1))).^a;
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF9(x,t)
    f=[];
    N=1+floor(10*abs(sin(0.5*pi*t)));
    g=ones(size(x,1),1);
    for i=2:size(x,2)
        tmp=x(:,i)-cos(4*t+x(:,1)+x(:,i-1));
        g=g+tmp.^2;
    end
    f1=g.*(x(:,1)+max(0, (0.1+0.5/N).*sin(2*N*pi*x(:,1))));
    f2=g.*(1-x(:,1)+max(0, (0.1+0.5/N).*sin(2*N*pi*x(:,1))));
    f=[f1,f2];
    xf=[x,f];
end

function xf= DF10(x,t)
    f=[];
    G=sin(0.5*pi*t);
    H=2.25+2*cos(0.5*pi*t);
    tmp=sin(2*pi*(x(:,1)+x(:,2)))/(1+abs(G));
    g=1+sum((x(:,3:end)-tmp).^2,2);
    f1=g.*sin(0.5*pi*x(:,1)).^H;
    f2=g.*sin(0.5*pi*x(:,2)).^H.*cos(0.5*pi*x(:,1)).^H;
    f3=g.*cos(0.5*pi*x(:,2)).^H.*cos(0.5*pi*x(:,1)).^H;
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF11(x,t)
    f=[];
    G=abs(sin(0.5*pi*t));
    g=1+G+sum((x(:,3:end)-0.5*G*x(:,1)).^2,2);
    y=pi*G/6+(pi/2-pi*G/3)*x(:,1:2);
    f1=g.*sin(y(:,1));
    f2=g.*sin(y(:,2)).*cos(y(:,1));
    f3=g.*cos(y(:,2)).*cos(y(:,1));
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF12(x,t)
    f=[];
    k=10*sin(pi*t);
    tmp1=x(:,3:end)-sin(t*x(:,1));
    tmp2=sin(floor(k*(2*x(:,1:2)-1))*pi/2);
    g=1+sum(tmp1.^2,2)+prod(tmp2,2);
    f1=g.*cos(0.5*pi*x(:,2)).*cos(0.5*pi*x(:,1));
    f2=g.*sin(0.5*pi*x(:,2)).*cos(0.5*pi*x(:,1));
    f3=g.*sin(0.5*pi*x(:,1));
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF13(x,t)
    f=[];
    G=sin(0.5*pi*t);
    p=floor(6*G);
    g=1+sum((x(:,3:end)-G).^2,2);
    f1=g.*cos(0.5*pi*x(:,1)).^2;
    f2=g.*cos(0.5*pi*x(:,2)).^2;
    f3=g.*sin(0.5*pi*x(:,1)).^2+sin(0.5*pi*x(:,1)).*cos(p*pi*x(:,1)).^2 +...
        sin(0.5*pi*x(:,2)).^2+sin(0.5*pi*x(:,2)).*cos(p*pi*x(:,2)).^2;
    f=[f1,f2,f3];
    xf=[x,f];
end

function xf= DF14(x,t)
     f=[];
     G=sin(0.5*pi*t);
     g=1+sum((x(:,3:end)-G).^2,2);
     y=0.5+G*(x(:,1)-0.5);
     f1=g.*(1-y+0.05*sin(6*pi*y));
     f2=g.*(1-x(:,2)+0.05*sin(6*pi*x(:,2))).*(y+0.05*sin(6*pi*y));
     f3=g.*(x(:,2)+0.05*sin(6*pi*x(:,2))).*(y+0.05*sin(6*pi*y));
     f=[f1,f2,f3];
     xf=[x,f];
 end

function xf= F8(x,t)
     f=[];
     G=sin(0.5*pi*t);
     H=1.25+0.75*sin(pi*t);
     gg=((x(:,1)+x(:,2))/2).^H;
     g=sum((x(:,3:end)-gg-G).^2,2);

     f1=(1+g).*cos(0.5*pi*x(:,2)).*cos(0.5*pi*x(:,1));
     f2=(1+g).*cos(0.5*pi*x(:,2)).*sin(0.5*pi*x(:,1));
     f3=(1+g).*sin(0.5*pi*x(:,2));

     f=[f1,f2,f3];
     xf=[x,f];
 end

 function xf= FDA4(x,t)
     f=[];
     G=sin(0.5*pi*t);
     g=sum((x(:,3:n)-G).^2, 2);


     f1=(1+g).*cos(0.5*pi*x(:,1)).*cos(0.5*pi*x(:,2));
     f2=(1+g).*cos(0.5*pi*x(:,1)).*sin(0.5*pi*x(:,2));
     f3=(1+g).*sin(0.5*pi*x(:,1));

     f=[f1,f2,f3];
     xf=[x,f];
 end

 function xf= FDA5(x,t)
     f=[];
     G=abs(sin(0.5*pi*t));
     F=1+100*sin(0.5*pi*t)^4;

     temp=x(:,3:n);
     g=G+sum((temp-G).^2, 2);


     f1=(1+g).*cos(0.5*pi*(x(:,1).^F)).*cos(0.5*pi*(x(:,2).^F));
     f2=(1+g).*cos(0.5*pi*(x(:,1).^F)).*sin(0.5*pi*(x(:,2).^F));
     f3=(1+g).*sin(0.5*pi*(x(:,1).^F));

     f=[f1,f2,f3];
     xf=[x,f];
 end

 function xf= FDA1(x,t)
     f=[];
     G=sin(0.5*pi*t);
     g=1+sum((x(:, 2:end)-G).^2, 2);
     f1=x(:,1);
     f2=g.*(1-sqrt(x(:,1)./g));
     f=[f1,f2];
     xf=[x,f];
 end

function xf= FDA2(x,t)
    f=[];
    n=size(x,2);

    f1=x(:,1);
    H=2*sin(0.5*pi*(t-1));
    temp=x(:,2:n-7);
    g=1+sum(temp.^2, 2);

    Htemp=H*ones(1,9);

    ntemp=x(n-8:end);
    sybc=H+sum((ntemp-Htemp/4).^2, 2);
    arbit=2^sybc;

    f2=g*(1-(f1/g).^arbit);
    f=[f1,f2];
    xf=[x,f];
end

function xf= FDA3(x,t)
    f=[];
    n=size(x,2);

    G=abs(sin(0.5*pi*t));
    F=10^(2*sin(0.5*pi*t));

    f1=0.5*(x(1)^F+x(2)^F);

    Gtemp=G*ones(1,n-2);

    temp=x(3:n);
    g=(1+G)+sum((temp-Gtemp).^2, 2);

    f2=g*(1-(f1/g).^0.5);
    f=[f1,f2];
    xf=[x,f];
end

function xf= DMOP1(x,t)
    f=[];
    n=size(x,2);
    H=0.75*sin(0.5*pi*t)+1.25;
    g=1+9*sum((x(2:end)).^2, 2)/(n-1);

    f1=x(1);
    f2=g*(1-(x(1)/(g^H)));
    f=[f1,f2];
    xf=[x,f];
end

function xf= DMOP2(x,t)
    f=[];
    G=abs(sin(0.5*pi*t));
    H=0.75*sin(0.5*pi*t)+1.25;
    g=1+sum((x(2:end)-G).^2, 2);
    f1=x(1);
    f2=g*(1-(x(1)/g)^H);
    f=[f1,f2];
    xf=[x,f];
end

function xf= DMOP3(x,t)
    f=[];
    n=size(x,2);
    G=abs(sin(0.5*pi*t));
    rdd=randperm(n);
    r=rdd(1);
    %
    %         r=1+floor((n-1)*G);
    tmp=setdiff(1:n,r);
    g=1+sum((x(tmp)-G).^2, 2);
    f1=x(r);
    f2=g*(1-(x(r)/g)^0.5);
    f=[f1,f2];
    xf=[x,f];
end

