%Kiran Rao
%ME 2016 - Section B
%Computer Project 1

function RaoKiranCP1()
    format long
    xVals = linspace(0.32,1,100); %generate x values for plot
    figure(1)
    plot(xVals, (((xVals)-0.32).^(0.32)+(xVals)-1)) %evaluate and plot y values
    xlabel('Velocity of Shockwave')
    ylabel('Various System Parameters')
    title('Relationship between Shockwave Velocity and System Parameters')
    grid on %turn on a grid for the plot
    function[rootValue, finalEa, numIterations] = secant(f, firstGuess, secondGuess, maxEa)
        format long %to print 8 digits in the results
        currentEa = 9e9; %set arbitrary high value for error
        numIterations = 0; 
        zeroTrue = fzero(f,(firstGuess+secondGuess)/2); %evaluate true zero of inputted function over interval
        numIterations = 0;
        guessVec = [];
        while maxEa < currentEa %until the current error is less than the stated maximum error
            nextGuess = secondGuess - ((f(secondGuess)*(firstGuess-secondGuess))/(f(firstGuess)-f(secondGuess))); %compute the next zero value estimate
            guessVec = [guessVec nextGuess]; %compile all past and current zero value estimates
            firstGuess = secondGuess; %shift estimates to evaluate new estimate
            secondGuess = nextGuess; %shift estimates to evaluate new estimate
            if numIterations>0 
                currentEa = 100*((abs(guessVec(end) - guessVec(end-1)))/(guessVec(end))); %calculate approximate percent relative error
            end
            numIterations = numIterations + 1; %increase counter for number of iterations by 1  
        end
        rootValue = nextGuess; %root estimate is the final estimated value after loop has ended
        finalEa = currentEa; %approximate percent relative error is final calculated error after loop ends
    
    end

[rootValue, finalEa, numIterations] = secant(@(x) exp(-x)-x,0,1,1); %evaluate example 6.6 in textbook using secant function
[rootValue2, finalEa2, numIterations2] = secant(@(x) (x-0.32).^(0.32)+(x)-1,0.32,1,1e-8);%evaluate shockwave function using secant function with error less than 1e-8
zeroTrueShock = fzero(@(x) (x-0.32).^(0.32)+(x)-1,.5); %find true zero of shockwave function

sprintf('Example 6.6 root approximation: %d \nExample 6.6 approximate percent relative error: %d \nExample 6.6 number of iterations: %d \nShockwave function root approximation: %d \nShockwave function approximate percent relative error: %d \nShockwave function number of iterations: %d \nShockwave function true root: %d',rootValue,finalEa,numIterations,rootValue2,finalEa2,numIterations2,zeroTrueShock)
end
