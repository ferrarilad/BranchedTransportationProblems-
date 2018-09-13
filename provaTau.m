clear all;
close all;
a = [5;2;1;.2];
b = [1;2;5;8];

x = linspace(0,10,1000);

tau = @(x) min(a*x+b);
figure(1)
plot(x, tau(x));

tauStar = zeros(size(x));
for i = 1:numel(x)
    j = find( x(i) > a , 1 );
    if isempty(j)
        tauStar(i)=-inf;
    else
        if(j-1)
            tauStar(i) = -(b(j)-b(j-1))*(x(i)-a(j-1))/(a(j)-a(j-1))-b(j-1);
        else
            tauStar(i) = 0;
        end
    end
end


figure(2)
plot(x,tauStar);

figure(3); hold on;
for i = 1:numel(x)    
    plot(x,x(i)*x-tauStar(i),'-r'); axis([0,10,0,20])

end
plot(x, tau(x));
hold off;


y = linspace(min(-b),max(-b),100);
xx=0:0.1:10;
out_1 = tauStar_1(-(xx+b(1)).^2,a,b);

figure(4)
plot(xx,out_1)


figure(5);  hold on;
for i = 1:numel(y)    
    plot(x,out_1(i)*x+(xx(i)+b(1))^2,'-r'); axis([0,10,0,20])

end
plot(x, tau(x));
hold off;