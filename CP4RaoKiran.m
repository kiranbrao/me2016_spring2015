%Kiran Rao
%ME 2016 - Section B
%902891012
%Computer Project 4

function CP4RaoKiran
    function[dYdt] = dY(t,Y) %calculate the derivative of an inital velocity and position vector using the differential equations used to model bungee jumping
        y1 = Y(1); %value of initial velocity (m/s)
        y2 = Y(2); %value of initial position (m)
        m = 70; %mass of jumper (kg)
        g = 9.81; %acceleration due to gravity (m/s^2)
        Cd = 0.2; %drag coefficient (kg/m)
        k = 50; %Hooke's spring constant of the bungee cord (N/m)
        lambda = 15; %damping constant of the cord (N-s/m)
        L = 10; %unstretched length of bungee cord
            if y2 < L %when the bungee cord has not yet stretched
                f1 =  (m*g - Cd*sign(y1)*(y1).^2)./m; %differential equation to equal to the rate of change of velocity
                f2 = y1;
                dYdt = [f1; f2]; %vector of acceleration as f1 and velocity as f2; the derivative of Y
            else %when the bungee cord begins to stretch
                f1 =  (m*g - Cd*sign(y1)*(y1).^2 - k*(y2 - L) - lambda*y1)./m; %differential equation to equal to the rate of change of velocity
                f2 = y1;
                dYdt = [f1; f2]; %vector of acceleration as f1 and velocity as f2; the derivative of Y
            end
    end

    tspan = 0:0.1:50; %vector representing time span of bungee jump in seconds 
    options = odeset('RelTol',1e-7); %set relative tolerance to 1 * 10^-7
    Y0 = [0 0]; %initial conditions are v(t=0) = 0m/s and x(t=0) = 0m
    [T,Y] = ode45(@dY,tspan,Y0,options); %use ode45 function to solve the differential equations; ode45 calls the dY function in order to solve the second order equation; output T is the original time vector, and output Y is an array with two columns containing the velocities and positions of the bungee jumper at the times in vector T
    
    figure(1)
    plot(T,Y(:,1),'k') %plot ode45 calculated velocity of bungee jump over time
    hold on
    
    function[Y] = RK4sys(tspan,Y0,f) %computes  4th Order Runge-Kutta method to solve the second order differential equations presented in function dY using the timespan and the intial velocity and position of the bungee jumper; output Y is an array with two columns containing the velocities and positions of the bungee jumper at the times in vector tspan
        Y = [];
        h = tspan(2) - tspan(1); %determine step size h
        k1 = f(tspan(1),Y0); %formula for intermediate step k1 in 4th Order Runge-Kutta equation using inital Y0 conditions
        k2 = f(tspan(1)+(h./2),Y0+(k1'.*h/2)); %formula for intermediate step k2 in 4th Order Runge-Kutta equation using inital Y0 conditions
        k3 = f(tspan(1)+(h./2),Y0+(k2'.*h/2)); %formula for intermediate step k3 in 4th Order Runge-Kutta equation using inital Y0 conditions
        k4 = f(tspan(1)+h,Y0+(k3'.*h)); %formula for intermediate step k4 in 4th Order Runge-Kutta equation using inital Y0 conditions
        Y = [Y0 ; Y0 + (1/6).*((k1' + 2.*k2' + 2.*k3' + k4').*(h))]; %formula for one step of 4th Order Runge-Kutta method to compute the next y-value; uses k1, k2, k3 and k4, and the intial Y0 conditions
        for i = 2:length(tspan)-1 %loop through the remaining time values to determine subsequent values of position and velocity using the 4th Order Runge Kutta method and the dY function
            k1 = f(tspan(i),Y(i,:)); %formula for intermediate step k1 in 4th Order Runge-Kutta equation
            k2 = f(tspan(i)+(h./2),Y(i,:)+(k1'.*h/2)); %formula for intermediate step k2 in 4th Order Runge-Kutta equation
            k3 = f(tspan(i)+(h./2),Y(i,:)+(k2'.*h/2)); %formula for intermediate step k3 in 4th Order Runge-Kutta equation
            k4 = f(tspan(i)+h,Y(i,:)+(k3'.*h)); %formula for intermediate step k4 in 4th Order Runge-Kutta equation
            Y = [Y ;(Y(i,:) + (1/6).*((k1' + 2.*k2' + 2.*k3' + k4').*(h)))]; %formula for one step of 4th Order Runge-Kutta method to compute the next y-value; uses k1, k2, k3 and k4
        end
    end
    
    f = @dY; %sets function handle f equal to function dY in order to paramterize second order differential equation into a first order differential equation
    [Y2] = RK4sys(tspan,Y0,f); %evaluate the positions and velocities of the bungee jumper over the times in tspan using the 4th Order Runge-Kutta method to solve the parametrized differential equations
    
    figure(1)
    plot(tspan,Y2(:,1),'r') %%plot RK4sys calculated velocity of bungee jump over time
    legend('Velocity of Bungee Jumper using ode45','Velocity of Bungee Jumper using RK4sys','Location','northeast') %create legend for figure 1
    xlabel('Time (s)') %labels x-axis
    ylabel('Velocity of Bungee Jumper (m/s)') %labels y-axis
    title('Velocity of Bungee Jumper over time using ode45 and RK4sys') %gives plot title
    
    maxFallode45 = max(Y(:,2)); %finds maxiumum value for position (fall length) of bungee jumper using the positions calculated from the ode45 function
    maxFallRK4sys = max(Y2(:,2)); %finds maxiumum value for position (fall length) of bungee jumper using the positions calculated from the RK4sys function
    
    sprintf('The maximum fall length calculated from ode45 is %f meters.\nThe maximum fall lengh calculated from RK4sys is %f meters.',maxFallode45,maxFallRK4sys)   
end
