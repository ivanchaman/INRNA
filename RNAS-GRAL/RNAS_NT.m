function RNAS_NT()
    tic;	
    clc; clear all;
    format long e
    fid=fopen('1_cos(110x).txt','w');
%entradas del algoritmo
    tol=1E-10; 
    %n=[5 10 20 30 40 50 60 70 80 90 100 200 300 500 750 1000 3000];           
    n=[5 11 21 31 41 51 61 71 81 91 101 201 301 501 751 1001 3001];
for t=1:17   
%Intervalos para el ajuste de la funcion
    efe=1;
    a=0; b=pi; h=pi/n(t);
    x=a:h:b;            
    c=zeros(n(t)+1,n(t)+1);        
    for j=1:n(t)+1         
        for i=1:n(t)+1
             c(j,i)=cos((j-1)*x(i));
        end
    end        
    w=rand(n(t)+1,1);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%intervalos de integracion
    a1=0; b1=pi; h1=(b1-a1)/n(t);
%Funciones Prueba
%Cambio de Variable Y=(a1+(b1-a1)/pi*x)    
    %F=sin(5*(a1+(b1-a1)/pi*x));       
    %F=sin(10*(a1+(b1-a1)/pi*x));
    %F=sin(50*(a1+(b1-a1)/pi*x));
    %F=sin(100*(a1+(b1-a1)/pi*x));

    %F=(a1+(b1-a1)/pi*x).^3.*sin(5*(a1+(b1-a1)/pi*x));
    %F=(a1+(b1-a1)/pi*x).^3.*sin(10*(a1+(b1-a1)/pi*x));
    %F=(a1+(b1-a1)/pi*x).^3.*sin(50*(a1+(b1-a1)/pi*x));
    %F=(a1+(b1-a1)/pi*x).^3.*sin(100*(a1+(b1-a1)/pi*x));

    %F=sin(5*cot((a1+(b1-a1)/pi*x))).*sin(2*(a1+(b1-a1)/pi*x));
    %F=sin(10*cot((a1+(b1-a1)/pi*x))).*sin(2*(a1+(b1-a1)/pi*x));
    %F=sin(50*cot((a1+(b1-a1)/pi*x))).*sin(2*(a1+(b1-a1)/pi*x));
    %F=sin(100*cot((a1+(b1-a1)/pi*x))).*sin(2*(a1+(b1-a1)/pi*x));    
    
    %F=exp((a1+(b1-a1)/pi*x)).*sin(5*(a1+(b1-a1)/pi*x));
    %F=exp((a1+(b1-a1)/pi*x)).*sin(10*(a1+(b1-a1)/pi*x));
    %F=exp((a1+(b1-a1)/pi*x)).*sin(50*(a1+(b1-a1)/pi*x));
    %F=exp((a1+(b1-a1)/pi*x)).*sin(100*(a1+(b1-a1)/pi*x));

    %F=sin(5*(a1+(b1-a1)/pi*x)).*cos(5*(a1+(b1-a1)/pi*x));
    %F=sin(5*(a1+(b1-a1)/pi*x)).*cos(10*(a1+(b1-a1)/pi*x));
    %F=sin(5*(a1+(b1-a1)/pi*x)).*cos(50*(a1+(b1-a1)/pi*x));
    %F=sin(10*(a1+(b1-a1)/pi*x)).*cos(50*(a1+(b1-a1)/pi*x));
    %F=sin(10*(a1+(b1-a1)/pi*x)).*cos(100*(a1+(b1-a1)/pi*x));

    %F=((1-cos((a1+(b1-a1)/pi*x))).^5).*sin(5*(a1+(b1-a1)/pi*x));
    %F=((1-cos((a1+(b1-a1)/pi*x))).^10).*sin(10*(a1+(b1-a1)/pi*x));
    %F=((1-cos((a1+(b1-a1)/pi*x))).^50).*sin(50*(a1+(b1-a1)/pi*x));
    %F=((1-cos((a1+(b1-a1)/pi*x))).^100).*sin(100*(a1+(b1-a1)/pi*x));        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Funciones con Singularidad
    %F=(a1+(b1-a1)/pi*x).*tan((a1+(b1-a1)/pi*x));
    %F=cos(3*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=cos(11*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=cos(21*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=cos(201*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=sin(3*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=(a1+(b1-a1)/pi*x).*cot(10*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*cot(60*(a1+(b1-a1)/pi*x)).^2;    
    %F=(a1+(b1-a1)/pi*x).*cot(110*(a1+(b1-a1)/pi*x)).^2;    
    %F=(a1+(b1-a1)/pi*x).*sec(10*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*sec(60*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*sec(110*(a1+(b1-a1)/pi*x)).^2;    
    %F=1./cos(10*(a1+(b1-a1)/pi*x));    
    %F=1./cos(60*(a1+(b1-a1)/pi*x));    
    F=1./cos(100*(a1+(b1-a1)/pi*x));            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %Algoritmo de Newton Truncado
    Y=c*w;
    error=F'-Y;    
    gk=-c*error;    %Gradiente gk
    hk=c*c';        %Hessiano  
    pj=zeros(n(t)+1,1);
    %pj(n(t)+1,1)=0 ;   %p(0)=0 
    dk=pj;
    rj=-gk;         %r(0)=-g(k)    
    dj=rj;          %d(0)=r(0)            
    deltaj=sum(rj.^2);  %delta(0)=sqrt(sum(r(0).^2))^2
    epsi=1e-10;
    k=0;
    l=0; 
    while efe > tol
    %while k < 2    
    alfak=1 ; 
    w=w+alfak*dk;
    while true
        qj=hk*dj;
        e1=sum(rj.^2);
        d2=dj'*qj      ;                                      
        d1=epsi*deltaj ;     
        if d2<=d1       %dj'*qj < e*deltaj                
            if l==0             %j=0
                dk=-gk;
            else                %j>0
                dk=pj;
            end
            break;
        end    
        alfaj=e1/d2 ;    %alfaj=sqrt(sum(rj.^2))^2/dj'*qj 
        pj=pj+alfaj*dj; %pj+1=pj+alfaj*dj          
        rj1=rj;
        rj=rj-alfaj*qj;  %rj+1=rj-alfaj*qj                                           
        al1=sqrt(sum(rj.^2));
        al2=sqrt(sum(gk.^2)) ;
        %if j>n+1
        if al1<=epsi*al2 || j>n(t)+1   %sqrt(sum(rj.^2))<niu*sqrt(sum(gk.^2)) || j>max
            dk=pj;
            break;
        end          
        c1=sum(rj.^2);
        c2=sum(rj1.^2);      
        betaj=c1/c2 ; %betaj=sqrt(sum(rj+1.^2))^2/sqrt(sum(rj.^2))^2                      
        dj=rj+betaj*dj; %dj+1=rj+1+betaj*dj        
        l=l+1;
      end         
      Y=c*w;
      error=F'-Y;
      efe=.5*sum(error.^2);                                                    
      k=k+1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   %Funcion Real    
    z=a1:h1:b1;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Integral de la Funci�n 1
    I2=(b1-a1)*w(1);
         for j=2:n(t)+1
             I2=I2+w(j)/(j-1)*(sin((j-1)*b1)-sin((j-1)*a1));
         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
%Integral de la Funcion 2
    k1=(b1-a1);
          I=k1*w(1);
          for j=2:n(t)+1
              c1=k1/((j-1)*pi);
              c2=sin((j-1)*pi) ;             
              c3=c1*c2;
              c4=w(j)*c3;
              I=I+c4;
          end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%Funcion aproximada
    for i=1:n(t)+1     
        Yx(i)=w(1);
        Yxx(i)=w(1);
        for j=2:n(t)+1          
            Yx(i)=Yx(i)+w(j)*cos((j-1)*pi/k1*(z(i)-a1));
            Yxx(i)=Yx(i)+w(j)*cos((j-1)*z(i));
        end            
    end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Resultados obtenidos
    %k
    %efe
    %I    
    %I2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Integral Real de las funciones 
    %real= (1-cos(50))/5;    
    %real= (1-cos(100))/10;
    %real= (1-cos(500))/50;
    %real= (1-cos(1000))/100;
    
    %real= ((300/25)-(6/5^4))*sin(50)+((60/5^3)-(1000/5)*cos(50));
    %real= (3-(6/10000))*sin(100)+((6/100)-(100)*cos(100));
    %real= ((3/25)-(6/6250000))*sin(500)+((6/12500)-(100/50)*cos(500));
    %real= ((3/100)-(6/100000000))*sin(1000)+((6/100000-10)*cos(1000));

    %real=5/2*pi*exp(-5);
    %real=5*pi*exp(-10);
    %real=25*pi*exp(-50);
    %real=50*pi*exp(-100);
    
    %real=5*exp(pi)/26+5/26;
    %real=10/101-10*exp(pi)/101;
    %real=50/2501-50*exp(pi)/2501;
    %real=100/10001-100*exp(pi)/10001;
    
    %real=sin(50)^2/10;
    %real=((sin(75)^2)/15)-((sin(25)^2)/25);
    %real=sin(275)^2/55-sin(225)^2/45;
    %real=sin(300)^2/60-sin(200)^2/40;
    %real=sin(550)^2/110-sin(450)^2/90;   
    
    %real=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Integral Real de Funciones con Singularidad
    %real=-pi*log(2);
    %real=-pi;
    %real=-pi;
    %real=pi;
    %real=pi;
    
    %real=0;
    
    %m=10;
    %m=60;
    m=110;
    %real=-4/m*cot(m*4)+log(sin(m*4))/m^2+8+2/m*cot(m*2)-log(sin(m*2))/m^2-2;
    %real=pi/m*tan(m*pi)+1/m^2*log(cos(m*pi));
    real=log(tan(pi/4+m*pi/2))/m;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%Error en las aproximaciones
%Error 1    
    e_I=abs(real-I);
%Error 2
    e_I2=abs(real-I2);
    fprintf('%d,%d,%9.20e,%9.20e,%9.20e,%9.20e,%9.20e\n',n(t),k,efe,I,I2,e_I,e_I2);
    fprintf(fid,'%d,%d,%9.20e,%9.20e,%9.20e,%9.20e,%9.20e\n',n(t),k,efe,I,I2,e_I,e_I2);
%     hold on
%     grid on
%     plot(z,Fz,'g');       
%     plot(z,Yx,'b');
%     plot(z,Yxx,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    tiempo=toc;
    fprintf(fid,'tiempo: %f\n',tiempo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    clear c;
    clear efe;
    clear F;
end 