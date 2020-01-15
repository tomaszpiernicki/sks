function [x]=korelLFM(a,b,fs)
%Funkcja oblicza funkcjê korelacji w sygna³ów a i b.
%Funkcja jest znormalizowana wzglêdem czêstotliwoœci próbkowania.
%Korelacja obliczana  jest metod¹ sukcesywnego mno¿enia wektorów.
%Zmienia automatycznie liczby zespolone na sprzê¿one tworz¹c wektory
%kolumnowe a' lub b' ,co jest potrzebne w korelacji funkcji zespolonych
%D³ugoœæ wektora w jest równa liczbie d³u¿szego z wektorów a,b
%Maksimum funkcji korelacji wystêpuje w po³owie krótszej funkcji

N=size(a,2);
M=size(b,2);

if N<=M
    w=zeros(1,M+N);
    b=[zeros(1,N) b zeros(1,N)];
    for m=1:M+N    
        w(m)=(b(m:N+m-1)*a')/fs;    %Mno¿enie wektorów a i b
    end
 x=w(round(N/2):end);
 xs=size(x);
 disp(xs)
else
    w=zeros(1,N+M);
    a=[zeros(1,M) a zeros(1,M)];
    for m=1:N+M    
        w(m)=(a(m:M+m-1)*b')/fs;    %Mno¿enie wektorów a i b
    end
    x=w(round(M/2):end);
end
end

