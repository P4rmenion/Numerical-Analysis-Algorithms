%Creating handle for sin(x) function.
f = @(x) sin(x);

%Lower bound.
a = 0;

%Upper bound.
b = pi/2;

%Number of sub intervals the original interval is divided in.
n = 11;

[I,e] = Simpson(f,a,b,n);

fprintf("[Simpson Method]\n");
fprintf("\nI = %.8f",I);
fprintf("\ne = %.8f",e);
fprintf("\n\nI refers to the approximated definite integral of f.");
fprintf("\ne refers to the corresponding error of the approximation.");

[I,e] = Trapezoidal(f,a,b,n);

fprintf("[Trapezoid Method]\n");
fprintf("\nI = %.8f",I);
fprintf("\ne = %.8f",e);
fprintf("\n\nI refers to the approximated definite integral of f.");
fprintf("\ne refers to the corresponding error of the approximation.");

function [I,e] = Simpson(f,a,b,n)
    %Formula ingredients [disecting interval and calculating border values]
    h = (b-a)/n;
    q = f(a) + f(b);
   
    %Calculating the second multiplication operand for Simpson's Qn formula.
    for i=1:2:n-1
        q = q+4*f(a+i*h);
    end

    for i=2:2:n-2
        q = q+2*f(a+i*h);
    end
    
    %Calculating final approximated integral.
    I = h/3*q;

    %Calculating error for Simpson method.
    e = abs(1 - I);
end

function [I,e] = Trapezoidal(f,a,b,n)
    %Disecting interval.
    h = (b - a)/n;
    x = a:h:b;
    I = 0;
    
    %Calculating integral according to trapezoidal rule.
    for i=1:1:n
        I = I + (x(i+1) - x(i))*(f(x(i+1)) + f(x(i)))/2;
    end
    
    e = abs(1 - I);
end