function [w]=korel(a,b,fs)
%Funkcja oblicza wartoœæ bezwzglêdn¹ funkcji korelacji w sygna³ów a i b.
%Funkcja jest znormalizowana wzglêdem czêstotliwoœci próbkowania.
%Korelacja obliczana  jest metod¹ sukcesywnego mno¿enia wektorów.
%Zmienia automatycznie liczby zespolone na sprzê¿one tworz¹c wektory
%kolumnowe a' lub b' ,co jest potrzebne w korelacji funkcji zespolonych
%D³ugoœæ wektora w jest równa liczbie d³u¿szego z wektorów a,b
%Wektor w nie zawiera stanu nieustalonego z pocz¹tku i koñca d³u¿szego
%wektora a lub b.

N=size(a,2);
M=size(b,2);

if N<=M
    w=zeros(1,M);
    b=[b zeros(1,N)];
    for m=1:M    
        w(m)=(b(m:N+m-1)*a')/fs;    %Mno¿enie wektorów a i b
    end
 
else
    w=zeros(1,N);
    a=[a zeros(1,M)];
    for m=1:N    
        w(m)=(a(m:M+m-1)*b')/fs;    %Mno¿enie wektorów a i b
    end
    
end
end

