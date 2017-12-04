%Kiran Rao
%ME 2016 - Section B
%902891012
%Homework 4: Problem 1

function HW4Pb1RaoKiran
    function [F] = dYdx(x,Y) %convert second order ODE into a first order ODE using an equivalent expression to subsitute for d2y/dx2; inputs are indpendent variable x and vector Y of [dy/dx and y]; output is vector F(x,y), which is the derivative of vector Y
        y1 = Y(1); %value of dy/dx
        y2 = Y(2); %value of y(x)
        f1 = (y2 - x + 2*y1)./7; %converts second order differential equation to first order by finding an equivalent expression for dy1/dx
        f2 = y1;
        F = [f1 ; f2]; %derivative of vector Y
    end
   
    ybound0 = 5; %lower bound
    ybound20 = 8; %upper bound
    guess1 = 0; %guess dy/dx(0) = 0
    guess2 = 1; %guess dy/dx(0) = 1
    x = linspace(0,20,100); %generate vector of x values between bounds
    Yvec10 = [guess1 ; ybound0]; %first guessed intiial condition for first order ODE solver
    Yvec20 = [guess2 ; ybound0]; %second guessed intiial condition for first order ODE solver

    options = odeset('RelTol', 1e-4);
    [X1,Y1] = ode45(@dYdx,x,Yvec10,options); %solve first order ODE using first guessed intiial conditions
    solvedY1 = Y1(:,2); %extract values for y(x) from first order ODE solver

    [X2,Y2] = ode45(@dYdx,x,Yvec20,options); %solve first order ODE using second guessed intiial conditions
    solvedY2 = Y2(:,2); %extract values for y(x) from first order ODE solver
    
    Y0 =(ybound20-solvedY1(end))./(solvedY2(end)-solvedY1(end))*(guess2 - guess1)+guess1; %solve for true solution to initial condition Y0 using linear interpolation of the first two guessed intial conditions and the final values of their respective solved y(x) column vectors

    Yvec0 = [Y0;ybound0]; %create true solution and use as intial condition
    [X3,Y3] = ode45(@dYdx,x,Yvec0,options); %find true y(x) values using the true solution as an intial condition to the first order ODE solver
    solvedY3=Y3(:,2); %extract values for y(x( from first order ODE solver
    
    figure (1)
    subplot(2,2,1) %create three plots for two guessed intial conditions and one solved initial condition
    plot(X1,solvedY1,'r', [0 20], [ybound0 ybound20], 'ko') %plot solution of first guessed intial condition from lower bound to upper bound; include true values of the boundary conditions
    xlabel('x') %labels x-axis
    ylabel('y') %labels y-axis
    title('Solution to problem 27.4 using initial guess y1(0) = 0') %gives plot title
    hold on
    subplot(2,2,2)
    plot(X2,solvedY2,'g', [0 20], [ybound0 ybound20], 'ko') %plot solution of second guessed intial condition from lower bound to upper bound; include true values of the boundary conditions
    title('Solution to problem 27.4 using initial guess y1(0) = 1') %gives plot title
    xlabel('x') %labels x-axis
    ylabel('y') %labels y-axis
    subplot(2,2,3)
    plot(X3,solvedY3,'b', [0 20], [ybound0 ybound20], 'ko')  %plot true solution as an intial condition from lower bound to upper bound; include true values of the boundary conditions
    title(sprintf('Solution to problem 27.4 using initial guess y1(0) = %f',Y0)) %gives plot title
    xlabel('x') %labels x-axis
    ylabel('y') %labels y-axis
end   