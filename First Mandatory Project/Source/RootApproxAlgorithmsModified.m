%Author: Parmenion Charistos
%Date: December 2020

%Description: Code for the second exercise of the first mandatory project for the course of Numerical Analysis,
%             given by the Aristotle university of Thessaloniki. Implementation of the modified Bisection, Newton - Raphson
%             and Secant methods for root approximation of a given function f.

%E2

%-----Modify the previous methods and recalculate roots-----

e = 0.000005;         %Desired precision.
                      
Mod_Bisection(e,0,3);   
Mod_Newton_Raphson(e,1);
Mod_Secant(e,1,1.5,0);


%Function Definitions

%-----Modified Bisection Method-----

function root = Mod_Bisection(e,a,b)
    N = 0;                           %Iterations counter.
    current_e = intmax;
    solved = true;

    while current_e > e
       N = N + 1;
       m = a + (b - a).*rand;        %New formula for updating interval (select random number within previous interval)
       if f(m)*f(a) < 0
           b = m;
       else
           if f(m)*f(b) < 0
               a = m;
           else
               disp("[Modified Bisection Method] Not a valid interval [a,b]. Initial product f(a)*f(b) should be < 0.")
               solved = false;
               break;
           end
       end
       current_e = abs(b-a);
       root = m;
    end

    if solved
        fprintf("\n\n[Modified Bisection Method] - Result: \n\n")   
        fprintf("Iterations: %d\n",N)                
        fprintf("Root: %.5f\n",root)     
        fprintf("f(%.5f) = %.15f\n",root,f(root)) 
    end
end


%-----Modified Newton Raphson Method-----

function root = Mod_Newton_Raphson(e,x_prev)
    N = 0;                                      %Iterations counter.
    current_e = intmax;                         %Current Error (initialized as max)

    while current_e > e
        N = N+1;                                                                                    %Incrementing iterations.
        x_next = x_prev - f(x_prev)/fd(x_prev) - ((f(x_prev).^2).*fdd(x_prev))/(2*(fd(x_prev).^3)); %Calculating next term of sequence.
        root = x_next;                                                                              %Updating current approxiimation.
        current_e = abs(x_next - x_prev);                                                           %Updating error after iteration.
        x_prev = x_next;                                                                            %Setting new term as previous term for the following iteration.
    end    

    fprintf("\n\n[Modified Newton Raphson Method] - Result: \n\n")   
    fprintf("Iterations: %d\n",N)                
    fprintf("Root: %.5f\n",root)   
    fprintf("f(%.5f) = %.15f\n",root,f(root)) 
end

%-----Modified Secant Method-----
               
function root = Mod_Secant(e,x0,x1,x2)
    N = 0;                      %Iterations counter.
    current_e = intmax;         %Current Error (initialized as max)

    while current_e > e
        N = N+1;                                                                %Incrementing iterations.
        q = f(x0)/f(x1);
        r = f(x2)/f(x1);
        s = f(x2)/f(x0);
        x3 = x2 - (r*(r-q)*(x2 - x1) + (1-r)*s*(x2 - x0))/((q-1)*(r-1)*(s-1));  %Calculating next term of sequence.
        root = x3;                                                              %Updating current approxiimation.
        current_e = abs(x3 - x0);                                               %Updating error after iteration.                                                                                 
        x2 = x1;
        x1 = x0;
        x0 = x3;
    end    

    fprintf("\n\n[Modified Secant Method Result]: \n\n")   
    fprintf("Iterations: %d\n",N)                
    fprintf("Root: %.5f\n",root) 
    fprintf("f(%.5f) = %.15f\n",root,f(root)) 
end

%--f function--
function y = f(x)
    y = 94*cos(x).^3 - 24*cos(x) + 177*sin(x).^2 - 108*sin(x).^4 - 72*cos(x).^3 * sin(x).^2 - 65;
end

%--f' function--
function y = fd(x)
    y = (216*cos(x).^2 - 432*cos(x))*sin(x).^3 + (-144*cos(x).^4 - 282*cos(x).^2 +354*cos(x) + 24) * sin(x);
end

%--f'' function--
function y = fdd(x)
    y = (432 - 432*cos(x))*sin(x).^4 + (1224*cos(x).^3 - 1296*cos(x).^2 + 564*cos(x) - 354)*sin(x).^2 - 144*cos(x).^5 - 282*cos(x).^3 + 354*cos(x).^2 + 24*cos(x);
end