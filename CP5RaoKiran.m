%Kiran Rao
%ME 2016 - Section B
%902891012
%Computer Project 5

function CP5RaoKiran
    L = 0.035;  %length of middle ear (m)
    H = 0.001; %height of cochlear ducts (m)
    Ko = 1e10; %stiffness at x = 0 (N/m^3)
    Mo = 0.5; %mass (kg)
    delta = 0.05; %damping constant 
    rho = 1000; %fluid density (kg/m^3)
    alpha = 3e2; %constant parameter(m^-1)
    omega = 2*pi*6000; %radian frequency 
    Zbm = @(x) (1i.*omega.*Mo + (delta.*(Ko.*Mo).^0.5 .* exp(-0.5.*alpha.*x)) + (Ko./(1i.*omega)).*exp(-1.*alpha.*x)); %impedance evaluated as a function of x
    
    xVec = linspace(0,L,1002); %vector of x values that will be used to evaluate pressure in the middle ear
    PL = 0; %initial condition P(L) = 0
    n = 1000; %number of iterations for matrix A
    A = sparse(n,n);  %create sparse(empty) matrix 
    F = sparse(n,1); 
    deltaX = xVec(2) - xVec(1); %calculate the step size between entries of the x vector
    c = @(x) (-2 - ((deltaX)^2*2i*rho*omega)./(H*Zbm(x))); %expression to simplify the system of equations created using the centered difference formula for the second derivative
    
    for i = 1:n %to populate vectors A and F
        A(i,i) = c(xVec(i+1)); %sets diagonal equal to c     
        if i<n
            A(i,i+1) = 1; %sets values right after diagonal equal to 1
        end
        if i>1
            A(i,i-1) = 1; %sets values right before diagonal equal to 1
        end
        
        F(i) = 0; %sets column vector F equal to 0
        if i==1
            F(1) = -deltaX; %sets the first entry in column vector F equal to step size delta X
            A(i,i) = c(xVec(i+1)) + 1; %sets first entry in matrix A equal to c + 1 to satisfy using forward difference for the first derivative
        end
    end
    
    P=A\F;  %solve to find vector P
    Po = P(1) - deltaX; %calculates P(0) using expression found for P(0) by expressing the first derivative as a forward difference
    Pextended = [Po ; P ; PL]; %adds P(0) and P(L) to calculated vector P
    figure(1)
    semilogy(xVec,abs(Pextended),'b') %plot the absolute value of the calculated pressure against the length of the middle ear using matrices of finite difference method
    vbm = (-2.*Pextended')./(Zbm(xVec)); %calculate basilar membrane velocity using calculated pressure vector
    figure(2)
    semilogy(xVec,abs(vbm),'k') %plot the absolute value of basilar membrane velocity against length of middle ear using matrices of finite difference method
    
    function[dPdx] = dX(x,P) %convert second order ODE into a first order ODE using an equivalent expression to subsitute for d2P/dx2; inputs are indpendent variable x and vector Y of [dy/dx and y]; output is vector F(x,y), which is the derivative of vector Y
        p1 = P(1); %value of dP/dx
        p2 = P(2); %value of P(x)
        dp1dx = ((2i.*rho.*omega.*p2)./(H.*Zbm(x))); %converts second order differential equation to first order by finding an equivalent expression for dp1/dx
        dPdx = [dp1dx ; p1]; %derivative of vector P
    end

    function [res] = BC(dP0,L) %residual vector must be set to 0 in order to satisyf boundary conditions
        dPdx0 = 1; %boundary condition dP/dx(0) = 1
        pL = 0; %Pressure at the apical end of the ear P(L) = 0      
        res(1)=dP0(1)-dPdx0;  % to enforce dP/dx(0) = 1
        res(2)=L(2)-pL;  % to enforce P(L) = 0
    end
  
    solinit=bvpinit(xVec,zeros(2,1));  %need to initialize the solution 
    sol=bvp4c(@dX,@BC,solinit);
    X=deval(sol,xVec);  %evaluate the solution on the mesh
    P2=X(2,:);  %Pressue values are the 2nd row of the matrix X
    
    figure(1)
    hold on
    semilogy(xVec,abs(P2),'g')  %plot the absolute value of the calculated pressure against the length of the middle ear using the bvp4c method
    legend('Using finite difference matrices','Using bvp4c function','Location','northeast') %create legend for figure 1
    ylabel('Fluid pressure (Pa)') %labels y-axis
    xlabel('Distance from basal end into cochlea (m)') %labels x-axis
    title('Fluid pressure of chochlea over the distance of the cochlea') %gives plot title
    
    vbm2 = (-2.*P2)./(Zbm(xVec)); %calculate basilar membrane velocity using second calculated pressure vector
    figure(2)
    hold on
    semilogy(xVec,abs(vbm2),'c') %plot the absolute value of basilar membrane velocity against length of middle ear using matrices of finite difference method
    legend('Using finite difference matrices','Using bvp4c function','Location','northeast') %create legend for figure 2
    ylabel('Basilar membrane velocity (m/s)') %labels y-axis
    xlabel('Distance from basal end into cochlea (m)') %labels x-axis
    title('Basilar membrane velocity over the distance of the cochlea') %gives plot title
end    