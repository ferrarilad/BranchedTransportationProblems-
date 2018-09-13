function y = tauStar_1(x,a,b)
y = zeros(size(x));
for i = 1:numel(x)
    j = find( x(i) > -b , 1 );
    if isempty(j)
        y(i)=-inf;
    else
        if(j-1)
            y(i) = -(a(j)-a(j-1))*(x(i)+b(j-1))/(b(j)-b(j-1))+a(j-1);
        else
            y(i) = a(1);
        end
    end
end

