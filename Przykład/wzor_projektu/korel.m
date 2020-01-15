function [w]=korel(a,b,fs)
%Funkcja oblicza warto�� bezwzgl�dn� funkcji korelacji w sygna��w a i b.
%Funkcja jest znormalizowana wzgl�dem cz�stotliwo�ci pr�bkowania.
%Korelacja obliczana  jest metod� sukcesywnego mno�enia wektor�w.
%Zmienia automatycznie liczby zespolone na sprz�one tworz�c wektory
%kolumnowe a' lub b' ,co jest potrzebne w korelacji funkcji zespolonych
%D�ugo�� wektora w jest r�wna liczbie d�u�szego z wektor�w a,b
%Wektor w nie zawiera stanu nieustalonego z pocz�tku i ko�ca d�u�szego
%wektora a lub b.

N=size(a,2);
M=size(b,2);

if N<=M
    w=zeros(1,M);
    b=[b zeros(1,N)];
    for m=1:M    
        w(m)=(b(m:N+m-1)*a')/fs;    %Mno�enie wektor�w a i b
    end
 
else
    w=zeros(1,N);
    a=[a zeros(1,M)];
    for m=1:N    
        w(m)=(a(m:M+m-1)*b')/fs;    %Mno�enie wektor�w a i b
    end
    
end
end

