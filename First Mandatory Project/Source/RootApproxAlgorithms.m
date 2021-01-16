%Author: Parmenion Charistos
%Date: December 2020

%Description: Code for the first exercise of the first mandatory project for the course of Numerical Analysis,
%             given by the Aristotle university of Thessaloniki. Implementation of the Bisection, Newton - Raphson
%             and Secant methods for root approximation of a given function f.

%E1

%-----Plot the function given in the interval [-2,2]-----

x = -2:0.1:2;
plot(x,f(x));

%-----Calculate the roots of function f in the selected interval-----

e = 0.000005;         %Desired precision.

Bisection(e,-1,2);    %Root Approx via Bisection.
Newton_Raphson(e,2);  %Root Approx via Newton_Raphson.
Secant(e,-2,2);       %Root Approx via Secant.

%Function Definitions

%-----Bisection Method-----

function root = Bisection(e,a,b)
    N = ceil(log((b-a)/e)/log(2));   %Number of iterations required.
    solved = true;

    for i = 1:N
       m = (a+b)/2;
       if f(m)*f(a) < 0
           b = m;
       else
           if f(m)*f(b) < 0
               a = m;
           else
               disp("[Bisection Method] Not a valid interval [a,b]. Initial product f(a)*f(b) should be < 0.")
               solved = false;
               break;
           end
       end
       root = m;
    end

    %Print number of iterations.
    %Print root found.
    
    if solved
        fprintf("[Bisection Method] - Result: \n\n")   
        fprintf("Iterations: %d\n",N)                
        fprintf("Root: %.5f\n",root)                 
        fprintf("f(%.5f) = %.15f\n",root,f(root)) 
    end
end


%-----Newton Raphson Method-----

function root = Newton_Raphson(e,x_prev)
    N = 0;                                      %Iterations counter.
    current_e = intmax;                         %Current Error (initialized as max)

    while current_e > e
        N = N+1;                                %Incrementing iterations.
        x_next = x_prev - f(x_prev)/fd(x_prev); %Calculating next term of sequence.
        root = x_next;                          %Updating current approximation.
        current_e = abs(x_next - x_prev);       %Updating error after iteration.
        x_prev = x_next;                        %Setting new term as previous term for the following iteration.
    end    

    fprintf("\n\n[Newton Raphson Method] - Result: \n\n")   
    fprintf("Iterations: %d\n",N)                
    fprintf("Root: %.5f\n",root)         
    fprintf("f(%.5f) = %.15f\n",root,f(root)) 
end


%-----Secant Method-----
    
function root = Secant(e,x_2prev,x_prev)           
    N = 0;                      %Iterations counter.
    current_e = intmax;         %Current Error (initialized as max)

    while current_e > e
        N = N+1;                                                                    %Incrementing iterations.
        x_next = x_prev - (f(x_prev)*(x_prev - x_2prev))/(f(x_prev) - f(x_2prev));  %Calculating next term of sequence.
        root = x_next;                                                              %Updating current approxiimation.
        current_e = abs(x_next - x_prev);                                           %Updating error after iteration.
        x_2prev = x_prev;
        x_prev = x_next;                                                            %Setting new term as previous term for the following iteration.
    end    

    fprintf("\n\n[Secant Method] - Result: \n\n")   
    fprintf("Iterations: %d\n",N)                
    fprintf("Root: %.5f\n",root) 
    fprintf("f(%.5f) = %.15f\n",root,f(root)) 
end

%--f function--
function y = f(x)
    y = exp(sin(x).^3) + x.^6 - 2*x.^4 - x.^3 - 1;
end

%--f' function--
function y = fd(x)
    y = 3*exp(sin(x).^3)*cos(x)*sin(x).^2+6*x.^5-8*x.^3-3*x.^2;
end

%--f'' function--
function y = fdd(x)
    y = 9*exp(sin(x).^3)*cos(x).^2*sin(x).^4 - 3*exp(sin(x).^3)*sin(x).^3 + 6*exp(sin(x).^3)*cos(x).^2*sin(x) + 30*x.^4 - 24*x.^2 - 6*x;
end


