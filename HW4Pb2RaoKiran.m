%Kiran Rao
%ME 2016 - Section B
%902891012
%Homework 4: Problem 2

function HW4Pb2RaoKiran
    deltaX = 2; %step size of x between each iteration 
    Y0 = 5; %lower boundary condition
    Y20 = 8; %upper boundary condition
    n = 9; %number of iterations
    
    A = sparse(n,n);  %create sparse(empty) matrix 
    F = sparse(n,1);
    a = (7/deltaX^2) - (1/deltaX); %expression to simplify the system of equations created using the centered difference formulas for first and second derivatives
    b = -((14/deltaX^2) + 1); %expression to simplify the system of equations created using the centered difference formulas for first and second derivatives
    c = (7/deltaX^2) + (1/deltaX); %expression to simplify the system of equations created using the centered difference formulas for first and second derivatives
    
    for i = 1:n
        A(i,i) = b; %sets diagonal equal to b       
        if i<n
            A(i,i+1) = a; %sets values right after diagonal equal to a
        end
        if i>1
            A(i,i-1) = c; %sets values right before diagonal equal to c
        end
        
        F(i) = -deltaX*i; %sets column vector F equal to respective x values at each iteration
        if i==1
            F(1) = -deltaX*i - c*Y0; %changes first entry of F to include lower boundary condition
        elseif i==n
            F(n) = -deltaX*i - a*Y20; %changes last entry of F to include upper boundary condition
        end
    end
    
    Y=A\F;  %solve to find vector T
         
    Textended=[Y0;Y;Y20]; %add the 2 boundary conditions to the vector Y
    x=linspace(0,20,n+2); %generate vector of x values

    plot(x,Textended,'b-') %plot vector Y with boundary conditions against vector of x values
    hold on
    xlabel('x')
    ylabel('Y')
    title('Solution to problem 27.5 using finite difference method (delta X = 2)')    
end                     