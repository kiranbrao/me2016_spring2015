%Kiran Rao
%ME 2016 - Section B
%902891012
%Computer Project 2

function CP2RaoKiran 
    function[a0,a1,rSquared] = linearRegression(xVals,yVals) % this function creates a least squares regression line for a set of data; inputs are x and y coordinates of data points; outputs are a0 and a1, used as the slope and intercept of the regression line, and rSquared, the coefficient of determination
        xSum = sum(xVals); %sum x data points
        ySum = sum(yVals); %sum y data points
        xyVals = xVals.*yVals; %multiply each x coordinate with the corresponding y coordinate
        xySum = sum(xyVals); %sum products of x and y coordinates
        n = length(xVals); %number of data points/entries used
        xSquared = xVals.^2; %square x data points
        xSquaredSum = sum(xSquared); %sum squares of x data points
        yAverage = ySum/n; %find y-bar or average value of y data points
        xAverage = xSum/n; %find x-bar or average value of x data points
        St = sum((yVals - ((yAverage)*(ones(1,n)))).^2); %calculate St using y values and y-bar
        a1 = ((n*xySum) - (xSum*ySum))/((n*xSquaredSum)-(xSum)^2); %calculate a1, the slope of the LSRL
        a0 = yAverage - a1*xAverage; %calculate a0, the intercept of the LSRL
        regressionEstimate = a0*ones(1,n) + a1*xVals; %calculate the estimated y-values for the LSRL
        Sr = sum((yVals - regressionEstimate).^2); %calculate Sr, the summation of the squares of the error
        rSquared = (St-Sr)/St; %determine coefficient of determination using St and Sr
    end
    displacement = 1:10; %x values of displacement for collected data
    ySpring1 = [4.5 5.5 8 12.3 15.5 19 21.5 24 25.5 30.3]; %y values of measured force for spring 1
    ySpring2 = [2.5 10 20 60 90 100 150 200 300 400];% y values of measured force for spring 2
    [a0Spring1,a1Spring1,rSquaredSpring1] = linearRegression(displacement,ySpring1); %create LSRL for spring 1
    [a0Spring2,a1Spring2,rSquaredSpring2] = linearRegression(displacement,ySpring2); %create LSRL for spring 2
    xRegression = linspace(0,10,200); %x values to plot the LSRL's of spring 1 and spring 2
    yRegSpring1 = a0Spring1 + a1Spring1*xRegression; %y values to plot LSRL of spring 1
    yRegSpring2 = a0Spring2 + a1Spring2*xRegression; %y values to plot LSRL of spring 2
    
    figure(1) 
    plot(displacement,ySpring1,'bx') %plot experimental data points for spring 1
    hold on
    plot(xRegression,yRegSpring1,'k--') %plot LSRL for spring 1
    xlabel('x (m)') %label x axis
    ylabel('Force (N)') %label y axis
    title('Spring 1') %give plot title
    
    figure(2)
    plot(displacement,ySpring2,'bx') %plot experimental data points for spring 2
    hold on
    plot(xRegression,yRegSpring2,'k--') %plot LSRL for spring 2
    xlabel('x (m)') %label x axis
    ylabel('Force (N)') %label y axis
    title('Spring 2') %give plot title
    
    function[a1,rSquared] = regressionNoIntercept(xVals,yVals) %this function creates a regression line where the intercept is at the origin (equal to 0); inputs are x and y coordinates of data points; outputs are a1, used as the slope of the regression line, and rSquared, the coefficient of determination
        ySum = sum(yVals); %sum y data points
        xyVals = xVals.*yVals; %multiply each x coordinate with the corresponding y coordinate
        xySum = sum(xyVals); %sum products of x and y coordinates
        n = length(xVals); %number of data points/entries used
        xSquared = xVals.^2; %square x data points
        xSquaredSum = sum(xSquared); %sum squares of x data points
        yAverage = ySum/n; %find y-bar or average value of y data points
        St = sum((yVals - ((yAverage)*(ones(1,n)))).^2); %calculate St using y values and y-bar
        a1 = (xySum/xSquaredSum); %calculate a1, the slope of the regression line
        regressionEstimate = a1*xVals; %calculate the estimated y-values for the regression line
        Sr = sum((yVals - regressionEstimate).^2); %calculate Sr, the summation of the squares of the error
        rSquared = (St-Sr)/St; %determine coefficient of determination using St and Sr
    end

    [a1NewSpring1,rSquaredNewSpring1] = regressionNoIntercept(displacement,ySpring1); %create zero intercept regression line for spring 1
    [a1NewSpring2,rSquaredNewSpring2] = regressionNoIntercept(displacement,ySpring2); %create zero intercept regression line for spring 2
    yRegNewSpring1 = a1NewSpring1*xRegression; %y values to plot zero intercept regression line of spring 1
    yRegNewSpring2 = a1NewSpring2*xRegression; %y values to plot zero intercept regression line of spring 2
    
    figure(1)
    hold on
    plot(xRegression,yRegNewSpring1,'r') %plot zero intercept regression line for spring 1
      
    figure(2)
    hold on
    plot(xRegression,yRegNewSpring2,'r') %plot zero intercept regression line for spring 2
    
    function[a,p,rSquared] = leastSquaresPowerLaw(xVals,yVals) %this function creates a power regression line for a set of data using the linearRegression function; inputs are x and y coordinates of data points; outputs are a, the coefficient of the power regression line, p, the exponent of the power regression line, and rSquared, the coefficient of determination
        xVals2 = log(xVals); %linearize x values based off calculations
        yVals2 = log(yVals); %linearize y values based off calculations
        [a0,a1] = linearRegression(xVals2,yVals2); %determine a0 and a1 (slope and intercept) of linearized power function 
        n = length(xVals); %number of data points/entries used
        ySum = sum(yVals); %sum y data points
        yAverage = ySum/n; %find y-bar or average value of y data points
        p = a1; %exponent of power regression line is the slope of linearized power function
        a = exp(a0); %coefficient of power regression line is the exponentiated intercept of the linearized power function
        f = @(x) a*x.^p; %create function handle for power regression function 
        regressionEstimate = f(xVals); %calculate the estimated y-values for the power regression line
        Sr = sum((yVals - regressionEstimate).^2); %calculate Sr, the summation of the squares of the error
        St = sum((yVals - ((yAverage)*(ones(1,n)))).^2); %calculate St using y values and y-bar
        rSquared = (St-Sr)/St; %determine coefficient of determination using St and Sr
    end
    
    [a,p,rSquaredPower] = leastSquaresPowerLaw(displacement,ySpring2); %create zero intercept regression function for spring 2
    f = @(x) a*x.^p; %create function handle for power regression function 
    yPowerRegression = f(xRegression); %y values to plot power regression line of spring 2
      
    figure(2)
    hold on
    plot(xRegression,yPowerRegression,'g') %plot power regression line for spring 2
    legend('Experimental Data', sprintf('y=a_1x+a_0, a_1=%.4f, a_0=%.4f, r^2=%.4f',a1Spring2,a0Spring2,rSquaredSpring2), sprintf('y=a_1x, a_1=%.4f, r^2=%.4f',a1NewSpring2,rSquaredNewSpring2), sprintf('y=ax^p, a=%.4f, p=%.4f, r^2=%.4f',a,p,rSquaredPower),'Location','northwest') %create legend for spring 2 plot
    set(gca,'XTick',0:2:10) %change increments of x axis
    
    figure(1)
    hold on
    legend('Experimental Data', sprintf('y=a_1x+a_0, a_1=%.4f, a_0=%.4f, r^2=%.4f',a1Spring1,a0Spring1,rSquaredSpring1), sprintf('y=a_1x, a_1=%.4f, r^2=%.4f',a1NewSpring1,rSquaredNewSpring1),'Location','northwest') %create legend for spring 1 plot
    set(gca,'XTick',0:2:10) %change increments of x axis
end   