function [p,hiba] = myRemez(n,epsilon,a,b)
    m = n+2;

    vege = -1;
    iter = 0;
    
    xk = zeros(1,m);
    for i = 1:m
        xk(i) = a + ((i-1) / (m-1)) * (b - a);
    end

    while (vege == -1 && iter < 10)
        B = f(xk);

        A = ones(m,m);
        for i = 1:(m/2)
            A(2*i,1) = -1;
        end
        for i = 1:m
            for j = 3:m
                A(i,j) = xk(i)^(j-2);
            end
        end 

        x = linsolve(A,B');

        p = zeros(1,n+1);
        for i = 1:(n+1)
            p(i) = x(n+2-i+1);
        end    

        g = @(x)(-1 * abs(f(x) - polyval(p,x)));
        kszi = fminbnd(g,a,b);
        
        if abs((-1 * g(kszi)) - x(1)) < epsilon
            vege = 1;
        else
            h = @(x)(f(x) * polyval(p,x));
            s = sign(h(kszi));
            if kszi < xk(2)
                if sign(h(xk(2))) == s
                    xk(2) = kszi;
                else
                    xk(n+1) = kszi;
                end
            elseif kszi > xk(n+1)
                if sign(h(xk(n+1))) == s
                    xk(n+1) = kszi;
                else
                    xk(2) = kszi;
                end
            else
                for i = 1:n+1
                    if kszi < xk(i+1)
                        ind = i;
                        break
                    end
                end
                if sign(h(xk(ind))) == s
                    xk(ind) = kszi;
                else
                    xk(ind+1) = kszi;
                end
            end
        end
    iter = iter + 1;
    end
    
    [z,hiba] = fminbnd(g,a,b);
    hiba = abs(hiba);
    
    xf = linspace(a,b,100);
    yf = f(xf);
    pf = polyval(p,xf);
    plot(xf,yf);
    hold on
    plot(xf,pf,'color','red');
    axis equal
    
    yg = abs(g(xf));
    plot(xf,yg,'color','green');
    hold off
end

