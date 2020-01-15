%oblLFM
disp('SYMULACJA SYSTEMU Z MODULACJ� LFM')
%Program funkcjonuje z programami komLFM.m i komLFM.fig
%Realizuje pr�bkowanie zwyk�e i kwadraturowe, 
%Program liczy stop� b��d�w w K transmisjach, kt�rych liczb� nale�y wpisa�.
%Przy du�ej liczbie K program wykonuje 8K rysunk�w!
%Nale�y "zakomentowa�" (% lub ctr r) instrukcje dorycz�ce wszystkich
%rysunk�w, a koniecznie rusunk�w do sprawozdania

%NADAJNIK
disp('Wyraz w kodzie ASCII')
sx=double(x);
disp(sx);
Mz=size(sx,2);       %Liczba znak�w
disp('Liczba znak�w');
disp(Mz);
disp('Liczba bit�w w transmisji');
M=Mz*8';
disp(M);
disp('Wyraz w zapisie binarnym');
Sb=dec2bin(sx,8);
disp(Sb)

%Utworzenie ciagu 0,1
S=zeros(Mz,8);
for m=1:Mz
    for k=1:8
        S(m,k)=str2double(Sb(m,k));
    end
end
disp('Wyraz jako macierz liczb binarnych')
disp(S);

%Rysunek sygna�u w bitach
sb=zeros(1,M);
for m=1:Mz
    for k=1:8
        n=(m-1)*8+k;
    sb(n)=S(m,k);
    end
end
skalan=1:M;
%Rysunek w GUI
axes(handles.axes1)
stem(skalan,sb,'markersize',2)
set(gca,'fontsize',8)
title('Bity nadane')
xlabel('n')
ylabel('sb(n)')
axis([0 M+1 0 1.5])

%Rysunek do sprawozdania
figure
set(gcf,'color','white')
stem(skalan,sb,'markersize',2)
set(gca,'fontsize',12)
title('Bity nadane')
xlabel('n')
ylabel('sb(n)')
axis([0 M+1 0 1.5])

%Generacja sygna��w LFM
disp('Cz�stotliwo�� pr�bkowania');
fs=4*f;
disp(fs);
disp('Liczba pr�bek w impulsie bitu')
Ni=round(T*fs);
disp(Ni);
disp('Czas trwania transmisji [s]')
disp(T*M);
disp('Liczba pr�bek w transmisji')
N=round(T*M*fs);
disp(N);

%Parametry sygna�u LFM
f1=f-0.7*B;       %Cz�stotliwo�� �rodkowa bitu 1
f0=f+0.7*B;      %Cz�stotliwo�� �rodkowa bitu 0
Bd=f-1.2*B;     %Cz�stotliwo�� dolna
Bg=f+1.2*B;    %Cz�stotliwo�� dolna
Bc=2.4*B;       %Szeroko�� pasma odbiornika

%Sygna�y wzorcowe 
a1=2*pi*f1/fs;
a2=2*pi*f0/fs;
b =2*pi*B/(2*fs);
n = 0:Ni-1; 
s1=sin((a1-b+b*n/Ni).*n);      %Sygna� bitu 1
s0=sin((a2-b+b*n/Ni).*n);      %Sygna� bitu 0
S1=abs(fft(s1));       %Widmo sygna�u bitu 1
S0=abs(fft(s0));       %Widmo sygna�u bitu 0
    
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wykres  widma sygna��w LFM (bit�w)
if p==1   %Pr�bkowanie zwyk�e
    %Rysunek w GUI
skalaf=(10^-3)*n*fs/Ni;
axes(handles.axes2)
plot(skalaf,S1,'k')
hold on
plot(skalaf,S0,'b')
set(gca,'fontsize',8)
set(gcf,'color','white')
xlabel('f [kHz]')
title('Widmo bit�w LFM')
hold off

%Rysunek do sprawozdania
figure
set(gcf,'color','white')
plot(skalaf,S1,'k')
hold on
plot(skalaf,S0,'b')
set(gca,'fontsize',12)
set(gcf,'color','white')
xlabel('f [kHz]')
title('Widmo bit�w LFM')
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modulacja
s=[];
for m=1:Mz
    for k=1:8
        if S(m,k)==1
        sb=s1;
        else
        sb=s0;
        end
        s=[s sb];     %Generowany sygna� w�skopasmowy
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sygna�y wzorcowe dla pr�bkowania kwadraturowego
if p==0
d=4*floor(f/Bc);       %Okres pr�bkowania kwadraturowego
Nik=floor(Ni/d);       %Liczba pr�bek kwadraturowych w bicie
disp('Cz�stotliwo�� pr�kowania kwadraturowego');
Fs=fs/d;
disp(Fs);
%Pr�bkowanie kwadraturowe 
for n=1:Nik
s1k(n)=s1((n-1)*d+1); %Sygna� bitu 1
s0k(n)=s0((n-1)*d+1); %Sygna� bitu 0
c1k(n)=s1((n-1)*d+2); %Sygna� bitu 1
c0k(n)=s0((n-1)*d+2); %Sygna� bitu 0
end
sk1=c1k+1i*s1k;       %Wzorzec bitu 1
sk0=c0k+1i*s0k;       %Wzorzec bitu 0
Sk1=abs(fft(sk1));     %Widmo sygna�u bitu 1
Sk0=abs(fft(sk0));     %Widmo sygna�u bitu 0

%Rysunek widma w GUI
skalafk=(10^-3)*(-Nik/2:Nik/2-1)*Fs/Nik;
axes(handles.axes2)
plot(skalafk,Sk1,'k')
hold on
plot(skalafk,Sk0,'b')
set(gca,'fontsize',8)
set(gcf,'color','white')
xlabel('f [kHz]')
title('Widmo bit�w LFM')
hold off

%Rysunek do sprawozdania
figure
set(gcf,'color','white')
plot(skalafk,Sk1,'k')
hold on
plot(skalafk,Sk0,'b')
set(gca,'fontsize',12)
set(gcf,'color','white')
xlabel('f [kHz]')
title('Widmo bit�w LFM')
hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ODBIORNIK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pr�bkowanie zwyk�e
if p==1 
disp('Liczba transmisji do wyznaczania stopy b��d�w')
K=input('K = ');
%P�tla do obliczania stopy b��d�w
for k=1:K  

%Szum
vr=0.5*10^(-snr/10);
sigma=sqrt(0.5*vr*fs/Bc)
noise=sigma*randn(1,N);
sn=s+noise;

%Korelacja dla bitu 1
 y1=abs(korelLFM(sn,s1,fs));
 %Rysunek do GUI
 skalat=10^3*(0:N-1)/fs;
 axes(handles.axes3)
 plot(skalat,10^3*y1(1:N))
 set(gca,'fontsize',8)
 xlabel('t [ms]')
 title('Sygna� bitu 1 po korealacji')
 axis([0 10^3*N/fs 0 1.2*max(10^3*y1)])
 
%Korelacja dla bitu 0
 y0=abs(korelLFM(sn,s0,fs));
 
 %Rysunek do sprawozdania
 figure
 set(gcf,'color','white')
 subplot(2,1,1)
 plot(skalat,10^3*y1(1:N))
 set(gca,'fontsize',12)
 xlabel('t [ms]')
 title('Sygna� bitu 1 po korealacji')
 axis([0 10^3*N/fs 0 1.2*max(10^3*y1)])
 subplot(2,1,2)
 plot(skalat,10^3*y0(1:N))
 set(gca,'fontsize',12)
 xlabel('t [ms]')
 title('Sygna� bitu 0 po korealacji')
 axis([0 10^3*N/fs 0 1.2*max(10^3*y1)]) 
 
  %Detekcja progowa
  p1=max(y1)*0.3;   %Pr�g detekcji
  z1=y1>p1;             %Sygna� bitu 1
  z0=y0>p1;             %Sygna� bitu 0
   
 %Sumowanie pr�bek w bitach  
   for n=1:M
    c1(n)=sum(z1((n-1)*Ni+(1:Ni) ));
    c0(n)=sum(z0((n-1)*Ni+(1:Ni) ));
   end
 %Detekcja progowa
 prog1=0.5*max(c1);
 prog0=0.5*max(c0);
 w1=c1>prog1;            %Bity 1
 w0=c0>prog0;            %Bity 0

 %Rysunek w GUI
 axes(handles.axes4)
 stem(w1,'markersize',2)
 set(gca,'fontsize',8)
 xlabel('n')
 title('Bity odebrane')
 axis([0 M+1 0 1.5])
 
 %Rysunek do sprawozdania
 figure
 set(gcf,'color','white')
  stem(w1,'markersize',2)
 set(gca,'fontsize',12)
 xlabel('n')
 title('Bity odebrane')
 axis([0 M+1 0 1.5])
 
  %Kontrola zgodno�ci bit�w widoczna 
  %Command Window
  for m=1:M
      if w1(n)==1&&w0(n)==1
          disp('B��D')
      elseif w1(n)==0&&w0(n)==0
      disp('B��D')
      end
  end
  
  %Konwersja na wyraz
  er=0;
for n=1:Mz
    m=(n-1)*8;
    cn=2^7*w1(m+1)+2^6*w1(m+2)+2^5*w1(m+3)+2^4*w1(m+4)+2^3*w1(m+5)+ 2^2*w1(m+6)+2^1*w1(m+7)+w1(m+8);
    
 %Kontrola zgodno�ci tekstu nadanego i odebranego
    if cn~=sx(n)
        er=er+1;
    end
    wyraz(n)=char(cn);  % Tekst odebrany
end
disp('Tekst odebrany');
disp(wyraz)
set(handles.y,'string',wyraz)
err(k)=er;   
end 
disp('Liczba b��d�w w transmisji')
disp(err)
%Sytopa b��d�w
stopa=sum(err)/K;
set(handles.error,'string',stopa)
end
%Koniec wersji z pr�bkowaniem zwykym

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pr�bkowanie kwadraturowe
if p==0   
    disp('Liczba transmisji do wyznaczania stopy b��d�w')
K=input('K = ');

for k=1:K  
%Szum
vr=0.5*10^(-snr/10);
sigma=sqrt(0.5*vr*fs/Bc)
noise=sigma*randn(1,N);
sn=s+noise;

%Pr�bkowanie kwadraturowe
Nk=round(N/d) ;      %Liczba pr�bek kwadraturowych
for n=1:Nk
s1k(n)=sn((n-1)*d+1); %Sk�adowa kwadraturowa
c1k(n)=sn((n-1)*d+2); %Sk�adowa synfazowa
end
sk=c1k+1i*s1k;       %Sygna� zespolony

%Korelacja dla bitu 1
  y1=abs(korelLFM(sk,sk1,Fs));
   
   %Rysunek w GUI
 skalatk=10^3*(0:Nk-1)/Fs;
 axes(handles.axes3)
 plot(skalatk,10^3*y1(1:Nk))
 set(gca,'fontsize',8)
 xlabel('t [ms]')
 title('Sygna� bitu 1 po korealacji')
axis([0 10^3*Nk/Fs 0 1.2*max(10^3*y1)])
  
%Korelacja dla bitu 0
   y0=abs(korelLFM(sk,sk0,Fs));
   
   %Rysunek do sprawozdania
 figure
 set(gcf,'color','white')
 subplot(2,1,1)
 plot(skalatk,10^3*y1(1:Nk))
 set(gca,'fontsize',12)
 xlabel('t [ms]')
 title('Sygna� bitu 1 po korealacji')
 axis([0 10^3*N/fs 0 1.2*max(10^3*y1)])
 subplot(2,1,2)
 plot(skalatk,10^3*y0(1:Nk))
 set(gca,'fontsize',12)
 xlabel('t [ms]')
 title('Sygna� bitu 0 po korealacji')
 axis([0 10^3*N/fs 0 1.2*max(10^3*y1)]) 

  %Detekcja progowa 
   p1=max(y1)*0.3;
   
   z1=y1>p1;     %Sygna� bit�w 1
   z0=y0>p1;     %Sygna� bit�w 0
     
   %Sumowanie pr�bek w bitach
   for n=1:M
    c1(n)=sum(z1((n-1)*Nik+(1:Nik) ));
    c0(n)=sum(z0((n-1)*Nik+(1:Nik) ));
   end
   
 %Detekcja progowa  
 prog1=0.5*max(c1);
 prog0=0.5*max(c0);
 w1=c1>prog1;              %Bity 1
 w0=c0>prog0;                %Bity 0
 
 %Rysunek w GUI
  axes(handles.axes4)
   stem(w1,'markersize',2)
   set(gca,'fontsize',8)
   xlabel('n')
  title('Bity odebrane')
   axis([0 M+1 0 1.5])
   
   %Rysunek do sprawozdania
   figure
   set(gcf,'color','white')
   stem(w1,'markersize',2)
   set(gca,'fontsize',12)
   xlabel('n')
  title('Bity odebrane')
   axis([0 M+1 0 1.5])
   
  %Kontrola zgodno�ci bit�w
  for m=1:M
      if w1(n)==1&&w0(n)==1
          disp('B��D')
      elseif w1(n)==0&&w0(n)==0
        disp('B��D')
      end
  end
  
  %Konwersja na wyraz
  er=0;
for n=1:Mz
    m=(n-1)*8;
    cn=2^7*w1(m+1)+2^6*w1(m+2)+2^5*w1(m+3)+2^4*w1(m+4)+2^3*w1(m+5)+ 2^2*w1(m+6)+2^1*w1(m+7)+w1(m+8);
    
 %Zliczanie b��d�w transmisji 
    if cn~=sx(n)
        er=er+1;
    end
    wyraz(n)=char(cn);
end
disp('Tekst odebrany');
disp(wyraz)
set(handles.y,'string',wyraz)
err(k)=er;
end 
disp('Liczba b��d�w w transmisji')
disp(err)
stopa=sum(err)/K;
set(handles.error,'string',stopa)
end



