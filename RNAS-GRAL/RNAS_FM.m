function RNAS_FM()
    tic
    clc; clear all
    format long e
%entradas del algoritmo
    tol=1E-10;     
    %n=[5 10 20 30 40 50 60 70 80 90 100 200 300 500 750 1000 ];           
    n=[5 11 21 31 41 51 61 71 81 91 101 201 301 501 751 1001];
for t=1:16   
%Intervalos para el ajuste de la funcion
    efe=1;
    lamda=1;
    a=0; b=pi; h=pi/n(t);
    x=a:h:b;            
    c=zeros(n(t)+1,n(t)+1);        
    k=0;
    for j=1:n(t)+1         
        for i=1:n(t)+1
             c(j,i)=cos((j-1)*x(i));
        end
    end        
    w=rand(n(t)+1,1);    
    r1=(sum(c.^2));
    r=sqrt(sum(r1));
    eta=1.35/r^2;   
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Funciones con Singularidad
    F=(a1+(b1-a1)/pi*x).*tan((a1+(b1-a1)/pi*x));
    %F=cos(3*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=cos(11*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=cos(21*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=cos(201*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=sin(3*(a1+(b1-a1)/pi*x))./cos((a1+(b1-a1)/pi*x));
    %F=(a1+(b1-a1)/pi*x).*cot(5*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*cot(10*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*cot(50*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*cot(100*(a1+(b1-a1)/pi*x)).^2;    
    %F=(a1+(b1-a1)/pi*x).*sec(5*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*sec(10*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*sec(50*(a1+(b1-a1)/pi*x)).^2;
    %F=(a1+(b1-a1)/pi*x).*sec(100*(a1+(b1-a1)/pi*x)).^2;
    %F=1./cos((a1+(b1-a1)/pi*x));    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %Algoritmo de Factor Momentum
    while efe > tol
    %while k < 500
        Y=c*w;
        error=F'-Y;
        w= w + eta*(lamda+1)*c'*error;
        efe=sqrt(sum(error.^2));
        efe=.5*efe^2;
        k=k+1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Funcion Real    
    z=a1:h1:b1;
    %Fz=sin(5*z);    
    %Fz=sin(10*z);
    %Fz=sin(50*z);
    %Fz=sin(100*z);

    %Fz=z.^3.*sin(5*z);
    %Fz=z.^3.*sin(10*z);
    %Fz=z.^3.*sin(50*z);
    %Fz=z.^3.*sin(100*z);

    %Fz=sin(5*cot(z)).*sin(2*z);
    %Fz=sin(10*cot(z)).*sin(2*z);
    %Fz=sin(50*cot(z)).*sin(2*z);
    %Fz=sin(100*cot(z)).*sin(2*z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Funciones con Sigularidad
    Fz=z.*tan(z);
    %Fz=cos(3*z)./cos(z);
    %Fz=cos(11*z)./cos(z);
    %Fz=cos(21*z)./cos(z);
    %Fz=cos(201*z)./cos(z);
    %Fz=sin(3*z)./cos(z);
    %Fz=z.*cot(5*z).^2;
    %Fz=z.*cot(10*z).^2;
    %Fz=z.*cot(50*z).^2;
    %Fz=z.*cot(100*z).^2;    
    %Fz=z.*sec(5*z).^2;
    %Fz=z.*sec(10*z).^2;
    %Fz=z.*sec(50*z).^2;
    %Fz=z.*sec(100*z).^2;
    %Fz=1./cos(z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %Fz=exp(z).*sin(5*z);
    %Fz=exp(z).*sin(10*z);
    %Fz=exp(z).*sin(50*z);
    %Fz=exp(z).*sin(100*z);

    %Fz=sin(5*z).*cos(5*z);
    %Fz=sin(5*z).*cos(10*z);
    %Fz=sin(5*z).*cos(50*z);
    %Fz=sin(10*z).*cos(50*z);
    %Fz=sin(10*z).*cos(100*z);

    %Fz=((1-cos(z)).^5).*sin(5*z);
    %Fz=((1-cos(z)).^10).*sin(10*z);
    %Fz=((1-cos(z)).^50).*sin(50*z);
    %Fz=((1-cos(z)).^100).*sin(100*z);
    
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
Yx=zeros(n(t)+1); 
Yxx=zeros(n(t)+1);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Integral Real de Funciones con Singularidad
    real=-pi*log(2);
    %real=-pi
    %real=-pi
    %real=pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
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
    
%Error en las aproximaciones
%Error 1    
    e_I=abs(real-I);
%Error 2
    e_I2=abs(real-I2);
    fprintf('%d,%d,%9.20e,%9.20e,%9.20e,%9.20e,%9.20e\n',n(t),k,efe,I,I2,e_I,e_I2);
    hold on
    grid on
    plot(z,Fz,'g');       
    plot(z,Yx,'b');
    plot(z,Yxx,'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    toc
    clear c;
    clear efe;
    clear F;
end 