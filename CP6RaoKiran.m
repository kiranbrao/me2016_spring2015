%Kiran Rao
%ME 2016 - Section B
%902891012
%Computer Project 6

function CP6RaoKiran
    R = 0.61803; %value of the Golden ratio R
    format long
    function[xopt,fxopt] = GoldenSectionSearch(f,xL,xR,Ea) %implements the Golden Section Search optimization method to find the maximum of a function between two bounds; inputs include a function handle, a left bound, a right bound, and the maxiumum allowable approximate relative percent error
        Eacurrent = 9e9; %set initial approximate relative percent error very high to initiate while loop
        while Eacurrent > Ea %while calculated approximate relative percent error is greater than the maxiumum allowable approximate relative percent error
            d = R *(xR - xL); %calculate d using left bound, right bound, and R
            x1 = xL + d; %calculate x1 using left bound and d
            x2 = xR - d; %calculate x2 using right bound and d
            f1 = f(x1); %evaluate x1 using function handle
            f2 = f(x2); %evaluate x2 using function handle
            if f2 > f1 %if the function evaluated at f2 is greater than the function evaluated at f1
               xR = x1; %x1 becomes the right bound for the next iteration
               xL = xL; %left bound remains the same
               xopt = x2; %the optimal/maximum value of x is considered to be x2 for this iteration
            else %if the function evaluated at f1 is greater than the function evaluated at f2
                xL = x2; %x2 becomes the left bound for the next iteration
                xR = xR; %right bound remains the same
                xopt = x1; %the optimal/maximum value of x is considered to be x1 for this iteration
            end 
            Eacurrent = 100 * (1-R) * abs((xL - xR)/xopt);
        end
        fxopt = f(xopt); %evaluate the optimal/maximum value of x using the function handle
    end

    V = 150000; %total number of units sold per year
    K1 = 10000; %(run)^-1
    K2 = 50; %(units-run)^1/2
    K3 = 300; %(run/unit)^2/3
    phi = @(S) ((V*(K1 + K2*S^(1/2)))/S) + K3*S^(2/3); %equation to model total cost of production and inventory
    negativePhi = @(S) -((V*(K1 + K2*S^(1/2)))/S) - K3*S^(2/3); %the above equation multiplied through by -1
    
    [SMin , costMinNeg] = GoldenSectionSearch(negativePhi,1,50000,1e-4); %use Golden Section Search function to find the maximum of the negative phi function, which will also be the minimum of the positive phi function
    costMin = -costMinNeg; %The minumum value of cost must be multiplied by -1 to make it positive again
    [SMin2, costMin2] = fminbnd(phi,1,50000); %use fminbnd function to also minimize the phi function and find optimal run size and minimum total cost; inputs include the function handle phi, a left bound, and a right bound
    
    % Using Golden Section Search, the optimal run size (S) is 15705.66083046405 units and the minimum value for total cost is $343497.4260153743
    % Using fminbnd, the optimal run size (S) is 15705.66142771371 units and the minimum value for total cost is $343497.4260153741
end
