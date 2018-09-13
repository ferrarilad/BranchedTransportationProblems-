function [f,df] = f_function(x,a,b) 
% a e b devono essere fatti in maniera che 
%     a(1) > a(2) > ... > a(N) >= 0 
% 0 = b(1) < b(2) < ... < b(N) 
% f = min( a(i)*x + b(i) )

f = a(1)*ones(size(x));
df =  zeros(size(x));
    for j = 1:length(x)
        i = find( -(x(j)-b(1))^2 >= -b , 1 );
        if i-1 
            f(j)  = -(a(i)-a(i-1))/(b(i)-b(i-1))*(-(x(j)-b(1)).^2+b(i-1))+a(i-1);
            df(j) = 2*(a(i)-a(i-1))/(b(i)-b(i-1))*(x(j)-b(1));
        end
    end
    
end