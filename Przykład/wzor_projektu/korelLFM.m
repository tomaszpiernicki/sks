function [x]=korelLFM(a,b,fs)
%Funkcja oblicza funkcj� korelacji w sygna��w a i b.
%Funkcja jest znormalizowana wzgl�dem cz�stotliwo�ci pr�bkowania.
%Korelacja obliczana  jest metod� sukcesywnego mno�enia wektor�w.
%Zmienia automatycznie liczby zespolone na sprz�one tworz�c wektory
%kolumnowe a' lub b' ,co jest potrzebne w korelacji funkcji zespolonych
%D�ugo�� wektora w jest r�wna liczbie d�u�szego z wektor�w a,b
%Maksimum funkcji korelacji wyst�puje w po�owie kr�tszej funkcji

N=size(a,2);
M=size(b,2);

if N<=M
    w=zeros(1,M+N);
    b=[zeros(1,N) b zeros(1,N)];
    for m=1:M+N    
        w(m)=(b(m:N+m-1)*a')/fs;    %Mno�enie wektor�w a i b
    end
 x=w(round(N/2):end);
 xs=size(x);
 disp(xs)
else
    w=zeros(1,N+M);
    a=[zeros(1,M) a zeros(1,M)];
    for m=1:N+M    
        w(m)=(a(m:M+m-1)*b')/fs;    %Mno�enie wektor�w a i b
    end
    x=w(round(M/2):end);
end
end

