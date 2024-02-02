% SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)


clc
clear all

% Define geometry
nDivX = 10 ;
nDivY = 10 ;

Lx = 1 ;
Ly = 1 ;

dx = Lx/nDivX ;
dy = Ly/nDivY ;

% Fluid properties
rho = 1000 ;
Re = 100 ;
U = 1 ;
L = 1 ;
mu = rho*U*L/Re ;

% Coordinates of the nodes
nCrds = (nDivX+1)*(nDivY+1);
x = linspace(0,Lx,nDivX+1);
y = linspace(0,Ly,nDivY+1);
[X,Y] = meshgrid(x,y);
crd = [reshape(X',nCrds,1) reshape(Y',nCrds,1)];

nPntX = nDivX+1 ;
nPntY = nDivY+1 ;

% Variable points
pres = ones(nPntX+1,nPntY+1);
U = zeros(nPntX, nPntY+1) ;
V = zeros(nPntX+1, nPntY) ;
pStar = 0*pres ;
UStar = 0*U ;
VStar = 0*V ;
pDash = 0*pres ;

% Time discretization
dt = 0.01 ;
iterMax = 500 ;
nSteps = iterMax ;

for step = 1:nSteps
    % Predict p*
    pStar = pres ;
    UStarPrev = UStar ;
    VStarPrev = VStar ;
    
    % Solve for u*
    for i=2:nPntX-1
        for j=2:nPntY
            vBar = 0.5*( VStarPrev(i,j) + VStarPrev(i+1,j) ) ;
            vDash = 0.5*( VStarPrev(i,j-1) + VStarPrev(i+1,j-1) );
            AStar  = -( (UStarPrev(i+1,j)^2 - UStarPrev(i-1,j)^2)/(2*dx) ...
                        + (UStarPrev(i,j+1)*vBar - UStarPrev(i,j-1)*vDash)/(2*dy) ) ...
                      + (mu/rho)*( (UStarPrev(i+1,j) - 2*UStarPrev(i,j) + UStarPrev(i-1,j))/(dx^2) + ...
                                   (UStarPrev(i,j+1) - 2*UStarPrev(i,j) + UStarPrev(i,j-1))/(dy^2));
            UStar(i,j) = UStarPrev(i,j) + AStar*dt - (dt/(rho*dx))*(pStar(i+1,j) - pStar(i,j)) ;
        end
    end
    % Satisfy boundary condition for UStar
    for i=1:nPntX
        UStar(i,1) = -UStar(i,2) ;         % U Bottom boundary 
        UStar(i,end) = 2-UStar(i,end-1);    % U Top boundary 
    end
    for j=1:nPntY+1
        UStar(1,j) = 0 ;                    % U Left boundary           
        UStar(end,j) = 0 ;                  % U Right boundary
    end
    UStar_inverted = flipud(UStar');
    
    % Solve for v*
    for i=2:nPntX
        for j=2:nPntY-1
            uDash = 0.5*( UStarPrev(i-1,j) + UStarPrev(i-1,j+1) ) ;
            uBar = 0.5*( UStarPrev(i,j) + UStarPrev(i,j+1) );
            BStar  = -( (VStarPrev(i+1,j)*uBar - VStarPrev(i-1,j)*uDash)/(2*dx) ...
                        + (VStarPrev(i,j+1)^2 - VStarPrev(i,j-1)^2)/(2*dy) ) ...
                      + (mu/rho)*( (VStarPrev(i+1,j) - 2*VStarPrev(i,j) + VStarPrev(i-1,j))/(dx^2) + ...
                                   (VStarPrev(i,j+1) - 2*VStarPrev(i,j) + VStarPrev(i,j-1))/(dy^2));
            VStar(i,j) = VStarPrev(i,j) + BStar*dt - (dt/(rho*dy))*(pStar(i,j+1) - pStar(i,j)) ;
        end
    end
    % Satisfy boundary condition for VStar
    for j=1:nPntY
        VStar(1,j) = -VStar(2,j) ;          % V Left boundary 
        VStar(end,j) = -VStar(end-1,j) ;    % V Right boundary 
    end
    for i=2:nPntX
        VStar(i,1) = 0 ;                    % V Bottom boundary       
        VStar(i,end) = 0 ;                  % V Top boundary
    end
    VStar_inverted = flipud(VStar');
    
    % Solve for pDash (Pressure Poisson equation)
    pDash = zeros(nPntX+1,nPntY+1);
    for iter=1:iterMax
        pDashOld = pDash ;
        for i=2:nPntX
            for j=2:nPntY
                a = 2*(dt/dx^2 + dt/dy^2) ;
                b = -dt/dx^2 ;
                c = -dt/dy^2 ;
                d = (rho/dx)*(UStar(i,j) - UStar(i-1,j)) + (rho/dy)*(VStar(i,j) - VStar(i,j-1)) ;
                %if (i~=2 || j~=2)
                pDash(i,j) = (-1/a) * (b*pDash(i+1,j) + b*pDash(i-1,j) + c*pDash(i,j+1) + c*pDash(i,j-1) + d) ;
                %end
                %cont(i,j) = d ;
            end
        end

        % Satisfy boundary condition for pressure
%     for j=2:nDivY+2
%         pDash(1,j) = pDash(3,j) ;     % P Left boundary
%         pDash(end,j) = pDash(end-2,j);% P Right boundary
%     end
%     for i=2:nDivX+2
%         pDash(i,1) = pDash(i,3);      % P Bottom boundary
%         pDash(i,end) = pDash(i,end-2);% P Top boundary
%     end
%     pDash(2,2) = 0;
    
        % Convergence criteria
        error_pDash = norm(pDash-pDashOld)./norm(pDash) ; 
        if (error_pDash < 1e-6)
            break;
        end 
    end
       
    % Update pressure
    alpha = 0.8 ;
    for i=2:nPntX
        for j=2:nPntY
            pres(i,j) = pStar(i,j) + alpha*pDash(i,j) ;
        end
    end        
    % Satisfy boundary condition for pressure
    for j=1:nPntY+1
        pres(1,j) = pres(2,j) ;     % P Left boundary
        pres(end,j) = pres(end-1,j);% P Right boundary
    end
    for i=1:nPntX+1
        pres(i,1) = pres(i,2);      % P Bottom boundary
        pres(i,end) = pres(i,end-1);% P Top boundary
    end
    %pres(2,2) = 0;
    
    % Convergence criteria
    errorP = norm(pDash); 
    errorU = norm(UStar - UStarPrev) ;
    errorV = norm(VStar - VStarPrev) ;
    fprintf('%d: ErrorU; ErrorV; ErrorP: %e %e %e\n',step,errorU, errorV, errorP);
    if (errorU < 1e-4 && errorV < 1e-4 && errorP < 1e-4)
        break;
    end
    
end

% Plot the solution
U = UStar ;
V = VStar ;
P = pres ;
for i=1:nPntX
    for j=1:nPntY
        uDisp(i,j) = 0.5*(U(i,j) + U(i,j+1)) ;
        vDisp(i,j) = 0.5*(V(i,j) + V(i+1,j)) ;
        pDisp(i,j) = 0.25*(P(i,j) + P(i+1,j) + P(i+1,j+1) + P(i,j+1)) ;
    end
end

xd = X;
yd = flipud(Y);

figure(1);
contourf(xd,yd,rot90(uDisp));
shading interp;
set(gca,'TickLabelInterpreter','latex','FontSize',30);
g1 = colorbar;
colormap('jet');
view([0 90]);
axis equal;
set(g1,'TickLabelInterpreter','latex','FontSize',30);
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u$','Interpreter','latex');

figure(2);
contourf(xd,yd,rot90(vDisp));
shading interp;
set(gca,'TickLabelInterpreter','latex','FontSize',30);
g1 = colorbar;
colormap('jet');
view([0 90]);
axis equal;
set(g1,'TickLabelInterpreter','latex','FontSize',30);
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$v$','Interpreter','latex');

figure(3);
contourf(xd,yd,rot90(pDisp));
shading interp;
set(gca,'TickLabelInterpreter','latex','FontSize',30);
g1 = colorbar;
colormap('jet');
view([0 90]);
axis equal;
set(g1,'TickLabelInterpreter','latex','FontSize',30);
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$p$','Interpreter','latex');

Re100_U = [1.0000	1.00000
0.9766	0.84123
0.9688	0.78871
0.9609	0.73722
0.9531	0.68717
0.8516	0.23151
0.7344	0.00332
0.6172	-0.13641
0.5000	-0.20581
0.4531	-0.21090
0.2813	-0.15662
0.1719	-0.10150
0.1016	-0.06434
0.0703	-0.04775
0.0625	-0.04192
0.0547	-0.03717
0.0000	0.00000];

%load('Re100_V.mat');

xDisp = [0:dx:Lx];
yDisp = [0:dy:Ly];

vDisp = rot90(vDisp);
uDisp = rot90(uDisp);

datav = vDisp(nDivX/2+1,:);
datau = uDisp(:,nDivY/2+1);

figure(4);
plot(xDisp,datav,'-r',Re100_V(:,1),Re100_V(:,2),'sb');

figure(5);
plot(flipud(datau),yDisp,'-r',Re100_U(:,2),Re100_U(:,1),'sb');


