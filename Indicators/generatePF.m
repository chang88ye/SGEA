function PF=generatePF(problem,t)

no=10000;  % expected number of PF points

switch problem
    case 1    %%%DF1
        x1=0:1/no:1;
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 2    %%%DF2
        x1=0:1/no:1;
        f1=x1;
        f2=1-x1.^0.5;
        PF=[f1;f2]';
    case 3    %%%DF3
        x1=0:1/no:1;
        H=1.5+sin(0.5*pi*t);
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 4   %%%DF4
        bt=1+abs(cos(0.5*pi*t));
        H=1.5+sin(0.5*pi*t);
        x1=0:1/no:bt;
        f1=x1.^H;
        f2=(bt-x1).^H;
        PF=[f1;f2]';
    case 5  %%%DF5
        x1=0:1/no:1;
        w=floor(10*sin(0.5*pi*t));
        f1=(x1+0.02*sin(w*pi*x1));
        f2=(1-x1+0.02*sin(w*pi*x1));
        PF=[f1;f2]';
    case 6  %%%DF6
        x1=0:1/no:1;
        at=0.2+2.8*abs(sin(0.5*pi*t));
        f1=(x1+0.1*sin(3*pi*x1)).^at;
        f2=(1-x1+0.1*sin(3*pi*x1)).^at;
        PF=[f1;f2]';
    case 7  %%%DF7
        x1=1:1/no:4;
        at=1+t;
        f1=at./x1;
        f2=x1./at;
        PF=[f1;f2]';
    case 8  %%%DF8
        x1=0:1/no:1;
        at=2.25+2*cos(2*pi*t);
        f1=(x1+0.1*sin(3*pi*x1));
        f2=(1-x1+0.1*sin(3*pi*x1)).^at;
        PF=[f1;f2]';
    case 9  %%%DF9 - disconnected
        x1=(0:1/no:1);
        nt=1+floor(10*abs(sin(0.5*pi*t)));
        f1=x1+max(0, (0.1+0.5/nt)*sin(2*nt*pi*x1));
        f2=1-x1+max(0, (0.1+0.5/nt)*sin(2*nt*pi*x1));
        k=f1+f2>1;
        f1(k)=[];
        f2(k)=[];
        PF=[f1;f2]';
    case 10 %%%DF10
        [x1,x2]=meshgrid(linspace(0,1,100));
        h=2.25+2*cos(0.5*pi*t);
        f1=sin(0.5*pi*x1).^h;
        f2=sin(0.5*pi*x2).^h.*cos(0.5*pi*x1).^h;
        f3=cos(0.5*pi*x2).^h.*cos(0.5*pi*x1).^h;
        PF=[f1(:),f2(:),f3(:)]; % column-wise concatenation by (:)
    case 11 %%%DF11
        [x1,x2]=meshgrid(linspace(0,1,100));
        gt=abs(sin(0.5*pi*t));
        y1=pi/6*gt+(pi/2-pi/3*gt)*x1;
        y2=pi/6*gt+(pi/2-pi/3*gt)*x2;
        f1=sin(y1);
        f2=sin(y2).*cos(y1);
        f3=cos(y2).*cos(y1);
        PF=[f1(:),f2(:),f3(:)];
    case 12 %%%DF12 -disconnected
        ObjNum=3;
        Dim=ObjNum;
        
        [x1,x2]=meshgrid(linspace(0,1,100));
        k=10*sin(pi*t);
        tmp2=abs(sin(floor(k*(2*x1-1))*pi/2).*sin(floor(k*(2*x2-1))*pi/2));
        g=1+tmp2;
        f1=g.*cos(0.5*pi*x2).*cos(0.5*pi*x1);
        f2=g.*sin(0.5*pi*x2).*cos(0.5*pi*x1);
        f3=g.*sin(0.5*pi*x1);

        P=[f1(:),f2(:),f3(:)];
        [PF,~]=NondomSort(P,0,ObjNum,Dim);
    case 13 %%%DF13 -disconnected
        ObjNum=3;
        Dim=ObjNum;
        [x1, x2]=meshgrid(linspace(0,1,100));
        G=sin(0.5*pi*t);
        p=floor(6*G);
        f1=cos(0.5*pi*x1).^2;
        f2=cos(0.5*pi*x2).^2;
        f3=(sin(0.5*pi*x1).^2+sin(0.5*pi*x1).*cos(p*pi*x1).^2 +...
            sin(0.5*pi*x2).^2+sin(0.5*pi*x2).*cos(p*pi*x2).^2);
        P=[f1(:),f2(:),f3(:)];
        [PF,~]=NondomSort(P,0,ObjNum,Dim);
    case 14 %%%DF14 -degenerate
        [x1,x2]=meshgrid(linspace(0,1,100));
        y1=0.5+sin(0.5*pi*t)*(x1-0.5);
        f1=(1-y1+0.05*sin(6*pi*y1));
        f2=(1-x2+0.05*sin(6*pi*x2)).*(y1+0.05*sin(6*pi*y1));
        f3=(x2+0.05*sin(6*pi*x2)).*(y1+0.05*sin(6*pi*y1));
        PF=[f1(:),f2(:),f3(:)];
    case 15 %%%F8
        [x1,x2]=meshgrid(linspace(0,1,100));
        f1=cos(0.5*pi*x1).*cos(0.5*pi*x2);
        f2=cos(0.5*pi*x1).*sin(0.5*pi*x2);
        f3=sin(0.5*pi*x1);
        PF=[f1(:),f2(:),f3(:)];
    case 16 %%%FAD4
        [x1,x2]=meshgrid(linspace(0,1,100));
        f1=cos(0.5*pi*x1).*cos(0.5*pi*x2);
        f2=cos(0.5*pi*x1).*sin(0.5*pi*x2);
        f3=sin(0.5*pi*x1);
        PF=[f1(:),f2(:),f3(:)];
    case 17 %%%FAD5
        [x1,x2]=meshgrid(linspace(0,1,100));
        gt=abs(sin(0.5*pi*t));
        ft=1+100*(sin(0.5*pi*t)^4);
        f1=(1+gt).*cos(0.5*pi*x1.^ft).*cos(0.5*pi*x2.^ft);
        f2=(1+gt).*cos(0.5*pi*x1.^ft).*sin(0.5*pi*x2.^ft);
        f3=(1+gt).*sin(0.5*pi*x1.^ft);
        PF=[f1(:),f2(:),f3(:)];
    case 18   %%%fad1
        x1=0:1/no:1;
        f1=x1;
        f2=1-x1.^0.5;
        PF=[f1;f2]';
    case 19  %%%fad2
        x1=0:1/no:1;
        h=1.5+sin(0.5*pi*t);
        H=2*h;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 20  %%%fad3
        x1=0:1/no:1;
        x2=0:1/no:1;
        gt=abs(sin(0.5*pi*t));
        ft=10^(2*sin(0.5*pi*t));
        f1=0.5*(x1.^ft+x2.^ft);
        f2=(1+gt).*(1-sqrt(f1./(1+gt)));
        PF=[f1;f2]';
    case 21   %%%dmop1
        x1=0:1/no:1;
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 22   %%%dmop2
        x1=0:1/no:1;
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x1;
        f2=1-x1.^H;
        PF=[f1;f2]';
    case 23   %%%dmop3
        x1=0:1/no:1;
        f1=x1;
        f2=1-x1.^0.5;
        PF=[f1;f2]';
    case 24   %%%F5
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';
    case 25   %%%F6
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';
    case 26   %%%F7
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';    
     case 27   %%%F9
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';   
     case 28   %%%F10
        x1=0:1/no:1;
        H=1.25+sin(0.5*pi*t)*0.75;
        f1=x1.^H;
        f2=(1-x1).^H;
        PF=[f1;f2]';   

    case 101,
        PF=[];
    
end
