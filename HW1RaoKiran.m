%Kiran Rao
%ME 2016 - Section B
%HW 1: Problem 4
%computation and plotting of true absolute error, true percent relative
%error, approximate absolute error, and approximate percent relative error
%for the Taylor series expansion of exp(x) for orders 0 through 15 where
%x=1

close all  %close all the figure windows
clear all  %clear all the variables
format long  % to print 8 digits in the results
currentApprox = 0;
trueAbsoluteError = []; %define empty vector for true absolute error
trueRelativeError = []; %define empty vector for true percent relative error
approxAbsoluteError = NaN; %the first term of approximate absolute error must be nonexistent because there is no n-1 for order 0
approxRelativeError = NaN; %the first term of approximate percent relative error must be nonexistent because there is no n-1 for order 0

for i=0:15 %orders 0 through 15 for taylor series approximation
    currentTerm = ((1)^i)./factorial(i); %taylor series approximation formula for exp(x) at x=1
    currentApprox = currentApprox + currentTerm; %adds previous terms to approximate nth term
    trueAbsoluteError = [trueAbsoluteError, (exp(1) - currentApprox)]; %create a vector of true absolute errors for orders 0 to 15
    trueRelativeError = [trueRelativeError, 100*((exp(1) - currentApprox)./exp(1))]; %create a vector of true percent relative errors for orders 0 to 15
    if i>0
        approxAbsoluteError = [approxAbsoluteError, (currentApprox - (currentApprox-currentTerm))]; %create a vector of approximate absolute errors for orders 1 to 15
        approxRelativeError = [approxRelativeError, 100*((currentApprox - (currentApprox-currentTerm))./currentApprox)]; %create a vector of approximate percent relative errors for orders 1 to 15
    end 
end

xOrders = 0:15; %create a vector of x values for orders 0 to 15

figure(1)
plot(xOrders,trueAbsoluteError) %plot true absolute error against orders 0 to 15 
xlabel('Order of Taylor Expansion')
ylabel('True Absolute Error')
title('True Absolute Error for Nth Order Expansion of exp(x)')
set(gca,'XTick',0:15) %change x-axis

figure(2)
plot(xOrders,trueRelativeError) %plot true percent relative error against orders 0 to 15 
xlabel('Order of Taylor Expansion')
ylabel('True Percent Relative Error')
title('True Percent Relative Error for Nth Order Expansion of exp(x)')
set(gca,'XTick',0:15) %change x-axis

figure(3)
plot(xOrders,approxAbsoluteError) %plot approximate absolute error against orders 1 to 15 
xlabel('Order of Taylor Expansion')
ylabel('Approximate Absolute Error')
title('Approximate Absolute Error for Nth Order Expansion of exp(x)')
set(gca,'XTick',0:15) %change x-axis

figure(4)
plot(xOrders,approxRelativeError) %plot approximate relative percent error against orders 1 to 15 
xlabel('Order of Taylor Expansion')
ylabel('Approximate Percent Relative Error')
title('Approximate Percent Relative Error for Nth Order Expansion of exp(x)')
set(gca,'XTick',0:15) %change x-axis

figure(5)
loglog(xOrders,approxRelativeError,'LineWidth',2) %plot approximate relative percent error against orders 1 to 15 on a logarithmic scale
xlabel('Order of Taylor Expansion')
ylabel('Approximate Percent Relative Error')
title('Approximate Percent Relative Error for Nth Order Expansion of exp(x) - Logarithmic')
set(gca,'XTick',[10^0 5 10^1 15 10^2]) %change x-axis
