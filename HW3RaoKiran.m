%Kiran Rao
%ME 2016 - Section B
%902891012
%Homework 3

function HW3RaoKiran
    format long %to print 8 digits in results
    function[ynew] = forwardEuler(t,y,h,f) %computes one step of forward Euler method to solve an initial value problem; t is the inital time, y is the inital y value, h is the step size, and f is the function used
        ynew = y + (f(t,y).*h); %formula for one step of forward Euler method to compute the next y-value
    end

    function[ynew] = RK4(t,y,h,f) %computes one step of 4th Order Runge-Kutta method to solve an initial value problem; t is the inital time, y is the inital y value, h is the step size, and f is the function used
        k1 = f(t,y); %formula for intermediate step k1 in 4th Order Runge-Kutta equation
        k2 = f(t+(h./2),y+(k1.*h/2)); %formula for intermediate step k3 in 4th Order Runge-Kutta equation
        k3 = f(t+(h./2),y+(k2.*h/2)); %formula for intermediate step k3 in 4th Order Runge-Kutta equation
        k4 = f(t+h,y+(k3.*h)); %formula for intermediate step k4 in 4th Order Runge-Kutta equation
        ynew = y + (1/6).*((k1 + 2.*k2 + 2.*k3 + k4).*(h)); %formula for one step of 4th Order Runge-Kutta method to compute the next y-value; uses k1, k2, k3 and k4
    end

    function[yV] = runODEsolver(f,y0,tV,method) %computes a vector of y-values at different times for an initial value problem given a function, an inital y-value, a vector of times, and a string specifying whether to use forward Euler method or 4th Order Runge-Kutta method
        h = tV(2) - tV(1); %calculates step size using two time values
        if strcmpi(method,'euler') 
            yV = [y0 forwardEuler(tV(1),y0,h,f)]; %compute first two entries of the vector of y-values with forward Euler method
            for i = 2:length(tV)-1 %index from the second time value to the end of the time vector
                yV = [yV forwardEuler(tV(i),yV(i),h,f)]; %compute remaining entries of vector of y-values with forward Euler method using remaining time intervals 
            end
        elseif strcmpi(method,'RK4')
            yV = [y0 RK4(tV(1),y0,h,f)]; %compute first two entries of the vector of y-values with 4th Order Runge-Kutta method
            for i = 2:length(tV)-1 %index from the second time value to the end of the time vector
                yV = [yV RK4(tV(i),yV(i),h,f)]; %compute remaining entries of vector of y-values with 4th Order Runge-Kutta method using remaining time intervals 
            end
        else
            sprintf('Last input must either be ''euler'' or ''RK4''') %prints instructions if incorrect string is used to specify either forward Euler method or 4th Order Runge-Kutta method
        end
    end

    f = @(t,y) 2*t - y^2; %create function handle for function in the intial value problem from homework problem 1
    y0 = 1; %initial y-value from initial value problem in homework problem 1
    tV = 0:0.1:10; %time interval to calculate y-values from initial value problem in homework problem 1 (0 seconds to 10 seconds with a step size of h = 0.1 seconds)
   
    yVEuler = runODEsolver(f,y0,tV,'euler'); %calculate y-values corresponding to specified time interval using forward Euler method and the intial value problem from homework problem 1
    yVRK4 = runODEsolver(f,y0,tV,'RK4'); %calculate y-values corresponding to specified time interval using 4th Order Runge-Kutta method and the intial value problem from homework problem 1
    
    figure(1)
    plot(tV,yVEuler,'b',tV,yVRK4,'r') %plot y-values calculated from the forward Euler method and 4th Order Runge-Kutta method and the initial value problem from homework problem 1
    xlabel('Time (s)') %labels x-axis
    ylabel('Value of y') %labels y-axis
    title('Value of y in Initial Value Problem over Time') %gives plot title
    legend('forward Euler method','4th Order Runge-Kutta method','Location','northwest') %create legend for figure 1
end