%Kiran Rao
%ME 2016 - Section B
%902891012
%Computer Project 3

function CP3RaoKiran 
    format long %to print 8 digits in the results
    function [I] = trapezoidal(xVals,yVals) %takes in vectors of x and y values then computes the area of trapezoids formed between each value to approximate the integral
        a = xVals(1); %left bound x-value
        b = xVals(end); %right bound x-value
        yVals2 = yVals; 
        yVals2(1) = []; 
        yVals2(end) = []; %get all internal y-values
        I = (b-a)*((yVals(1) + 2*sum(yVals2) + yVals(end))/(2*(length(yVals)-1))); %use composite trapezoidal formula to calculate integral approximation
    end

    Itrap = [];
    for n = 1:1000 %loop through cases where there are 1 to 1000 segements used
        xValsTrap = linspace(0,pi,n+1); %space out x-values according to number of segments
        yValsTrap = sin(xValsTrap); %calculate corresponding y-values
        Itrap = [Itrap trapezoidal(xValsTrap,yValsTrap)]; %vector of composite trapezoidal integral approximations for cases where the number of segments ranges from 1 to 1000
    end
    
    truePercentRelativeTrap = 100*(abs(2-Itrap)/2); %calculate true percent relative error for composite trapezoidal approximations
    xValsN = 1:1000; 

    figure(1)
    loglog(xValsN,truePercentRelativeTrap,'g') %plot true percent relative error against the number of segments in the composite trapezoidal integral approximation
    xlabel('Number of Segments') %labels x-axis
    ylabel('True Percent Relative Error') %labels y-axis
    title('True Percent Relative Error vs Number of Segments') %gives plot title
    hold on

    function[I] = compositeSimp13(xVals,yVals) %computes composite Simpson's 1/3 Rule for vectors of x-values and y-values using a quadratic interpolation
        a2 = xVals(1); %left bound x-value
        b2 = xVals(end); %right bound x-value
        yVals3 = yVals;
        yVals3(1) = [];
        yVals3(end) = []; %get all internal y-values
        I = (b2-a2)*((yVals(1) + 4*sum(yVals3(1:2:end)) + 2*sum(yVals3(2:2:end)) + yVals(end))/(3*(length(yVals)-1))); %use composite Simpson's 1/3 formula to calculate integral approximation
    end
    
    I13 = [];
    for n = 2:2:1000 %loop through cases where there are an even number of segments because Simpson's 1/3 rule only works with an even number of segments
        xVals13 = linspace(0,pi,n+1); %space out x-values according to number of segments
        yVals13 = sin(xVals13); %calculate corresponding y-values
        I13 = [I13 compositeSimp13(xVals13,yVals13)]; %vector of composite Simpson's 1/3 Rule integral approximations for cases where the number of segments ranges from 2 to 1000 in even intervals
    end
      
    truePercentRelative13 = 100*(abs(2-I13)/2); %calculate true percent relative error for composite Simpson's 1/3 approximations
    xValsN2 = 2:2:1000; 
    figure(1)
    loglog(xValsN2,truePercentRelative13,'rx') %plot true percent relative error against the number of segments in the composite Simpson's 1/3 Rule integral approximation

    function[I] = Simp38(h,fVals) %computes a single application Simpson's 3/8 Rule integral approximation for a specifed segment width and corresponding f(x) values
        I = (3*h/8)*(fVals(1) + 3*fVals(2) + 3*fVals(3) + fVals(4)); %formula for a single application of Simpson's 3/8 Rule
    end

    function[I] = SimpInt(xVals,yVals) %computes an approximated integral for vectors of x-values and y-values using a different Simpson's rule depending on the number of segments used
        if length(xVals) == 2 
            I = trapezoidal(xVals,yVals); %when there is only one segment, use trapezoidal rule  
        elseif any(length(xVals) == 1:2:1001)
            I = compositeSimp13(xVals,yVals); %when there is an even number of segments (the length of the x-vector is odd), use the composite Simpson's 1/3 Rule
        elseif length(xVals) == 4
            h = ((xVals(end)-xVals(1))/3); 
            fVals = yVals;
            I = Simp38(h,fVals); %when there are exactly 3 segments, use the single application Simpson's 3/8 Rule
        elseif any(length(xVals) == 6:2:1000)
            h = ((xVals(end)-xVals(end-3))/3);
            fVals = yVals(end-3:end);
            I = compositeSimp13(xVals(1:end-3),yVals(1:end-3)) + Simp38(h,fVals); %when there are an ood number of segments, the approximated integral should be a sum of the single application of Simpson's 3/8 Rule for the last 3 segments and the composite Simpson's 1/3 rule for the rest of the segments (an even number of segments)
        end
    end

    IAll = [];
    for n = 1:1000 %loop through cases where there are 1 to 1000 segements used
        xValsAll = linspace(0,pi,n+1); %space out x-values according to number of segments
        yValsAll = sin(xValsAll); %calculate corresponding y-values
        IAll = [IAll SimpInt(xValsAll,yValsAll)]; %vector of integral approximations using the SimpInt function for cases where the number of segments ranges from 1 to 1000
    end
    
    truePercentRelativeAll = 100*(abs(2-IAll)/2); %calculate true percent relative error for the Simpson's rules used in the SimpInt function
    figure(1)
    loglog(xValsN,truePercentRelativeAll,'b') %plot true percent relative error against the number of segments in the SimpInt function integral approximation 
    legend('trapezoidal','compositeSimp13', 'SimpInt','Location','northeast') %create legend for figure 1
    
    loading = importdata('loading.dat'); %import data from loading.dat file
    strain = loading(1:end,3); %index strain values from loading.dat
    stress = loading(1:end,2); %index stress values from loading.dat
    figure(2)
    plot(strain,stress) %plot stress vs strain from loading.dat
    xlabel('Strain') %labels x-axis
    ylabel('Stress (MPa)') %labels y-axis
    title('Stress-Strain Curve of Viscoelastic Material') %gives plot title
    
    workTrap = trapezoidal(strain,stress); %approximate the integral from time = 0 to time = 2 for the stress-strain curve using a composite trapezoidal calculation 
    workSimp = SimpInt(strain,stress); %approximate the integral from time = 0 to time = 2 for the stress-strain curve using the SimpInt function 
end