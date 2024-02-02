clc 
clear 

%% setting varibles 
tic
xgridsize = 21;
ygridsize = 21;
dx = 1/(xgridsize-1);
dy = 1/(ygridsize-1);
dt =  0.01;
rho = 1000;
mew = 10;
pfeild = ones(ygridsize + 1,xgridsize + 1);
pfinal = zeros(ygridsize ,xgridsize );
ufinal = zeros(ygridsize ,xgridsize );
vfinal = zeros(ygridsize ,xgridsize );
pcorr = zeros(ygridsize + 1,xgridsize + 1);
pcorrnew = zeros(ygridsize + 1,xgridsize + 1);
ufeild = zeros(ygridsize + 1 ,xgridsize);
vfeild = zeros(ygridsize,xgridsize  + 1);
unewfeild = zeros(ygridsize + 1,xgridsize);
vnewfeild = zeros(ygridsize,xgridsize + 1);
tollarance = 0.0001;
iterations = 0;
pnorm = 1;
dpnorm = 1;
parray = [];

%% boundary conditions repeat code 
% for i = 1:xgridsize
%     ufeild(1,i) = 2 - ufeild(1,i);
%     unewfeild(1,i) = 2 - unewfeild(2,i);
%     ufeild(ygridsize + 1,i) = - ufeild(ygridsize,i);
%     unewfeild(ygridsize + 1,i) = - unewfeild(ygridsize,i);
% end
% 
% for i = 1:ygridsize
%     vfeild(i,1) = - vfeild(i,2);
%     vfeild(i,xgridsize + 1) = - vfeild(i , xgridsize);
%     vnewfeild(i,1) = - vnewfeild(i,2);
%     vnewfeild(i,xgridsize + 1) = - vnewfeild(i , xgridsize);
% end
%% Non repeatable code 

for i = 1:ygridsize + 1
    ufeild(i,1) = 0;
    ufeild(i,xgridsize) = 0;
    unewfeild(i,1) = 0;
    unewfeild(i,xgridsize) = 0;
end

for i = 1:xgridsize + 1
    vfeild(1,i) = 0;
    vfeild(ygridsize,i) = 0;
    vnewfeild(1,i) = 0;
    vnewfeild(ygridsize,i) = 0;
end

%% non linear iterative loop
 for n2 = 1:100000
    iterations = iterations + 1
    %calculating u*
    for i = 2:ygridsize
        for j = 2 : xgridsize - 1 
        v1 = 0.5 * (vfeild(i-1,j) + vfeild(i-1,j+1));
        v2 = 0.5 * (vfeild(i,j) + vfeild(i,j+1));
        a = -( (((rho * ufeild(i,j+1) * ufeild(i,j+1)) - (rho * ufeild(i,j-1)*ufeild(i,j-1)))/(2 * dx))...
            + (((rho * ufeild(i-1,j)*v1) - (rho*ufeild(i+1,j)*v2)) / (2*dy)))...
            + mew*( ((ufeild(i,j+1)-2*ufeild(i,j)+ufeild(i,j-1))/(dx*dx)) + ...
            ((ufeild(i+1,j) - 2*ufeild(i,j) + ufeild(i-1,j))/(dy*dy)) );
        unewfeild(i,j) = ufeild(i,j) + (1/rho)*(a*dt -(dt/dx)*(pfeild(i,j+1)-pfeild(i,j)));
        end
        end
    
    %code to be repeated after every update of u
    for i = 2:xgridsize-1
    %ufeild(1,i) = 2 - ufeild(2,i);
    unewfeild(1,i) = 2 - unewfeild(2,i);
    %ufeild(ygridsize + 1,i) = - ufeild(ygridsize,i);
    unewfeild(ygridsize + 1,i) = - unewfeild(ygridsize,i);
    end

    %calculating v*

    for i = 2:ygridsize -1 
        for j = 2 : xgridsize
        u1 = 0.5 * (ufeild(i,j-1)+ufeild(i + 1,j-1)); % at point c
        u2 = 0.5 * (ufeild(i,j)+ufeild(i + 1,j)); % at point d 
        b = - ((( (rho*u2*vfeild(i,j+1)) - (rho*u1*vfeild(i,j-1)) )/(2*dx)) + (( (rho * vfeild(i-1,j)*vfeild(i-1,j)) - (rho * vfeild(i+1,j)*vfeild(i+1,j)))/(2*dy))) + mew*(((vfeild(i,j+1) - 2*vfeild(i,j) + vfeild(i,j-1))/(dx*dx))+((vfeild(i+1,j) - 2*vfeild(i,j) + vnewfeild(i-1,j))/(dy*dy)));
        vnewfeild(i,j) = vfeild(i,j) + (1/rho) * (b*dt - (dt/dy)*(pfeild(i,j) - pfeild(i+1,j)));
        end 
    end

    %code to be repeated after every update of v
    for i = 2:ygridsize-1
        %vfeild(i,1) = - vfeild(i,2);
        %vfeild(i,xgridsize + 1) = - vfeild(i , xgridsize);
        vnewfeild(i,1) = - vnewfeild(i,2);
        vnewfeild(i,xgridsize + 1) = - vnewfeild(i , xgridsize); 
    end 
    
    
    %solving for pressure correction values 
    pcorrnew = zeros(ygridsize + 1,xgridsize + 1);
    for n1 = 1:100000
        
        pcorr =pcorrnew ;
        for i = 2:ygridsize
            for j = 2:xgridsize
               a = 2 * (((dt)/(dx*dx)) + ((dt)/(dy*dy)));
               b = -((dt)/(dx*dx));
               c = -((dt)/(dy*dy));
               d = (1/dx) * ((rho * unewfeild(i,j)) - (rho*unewfeild(i,j-1))) + (1/dy)*((rho*vnewfeild(i-1,j))...
                   - (rho*vnewfeild(i,j)));
               pcorrnew(i,j) = (-1/a) * ( b * pcorrnew(i,j+1) + b * pcorrnew(i,j-1) + c * pcorrnew(i+1,j) + c*pcorrnew(i-1,j) + d);
%                if i == ygridsize - 1 && j == 2
%                    pcorrnew(i,j) = 0;
%                end
            end
        end
        
    dpnorm = norm(pcorrnew - pcorr) / norm(pcorrnew);

    if dpnorm <1e-6  && n1 > 2
        break
    end 
    end

    %end of loop 
    
    
    %update pressure values 
    alpha = 0.1;
    pfeild = pfeild + alpha.*pcorrnew;

    %setting newman boundary conditions
        for i = 1:xgridsize + 1
            pfeild(1,i) = pfeild(2,i);
            pfeild(ygridsize + 1 ,i) = pfeild(ygridsize,i);
        end 
        for i = 1:ygridsize + 1
            pfeild(i,1) = pfeild(i,2);
            pfeild(i,xgridsize + 1) = pfeild(i,xgridsize);
        end 

    pnorm = norm(pcorrnew);
    parray = [parray,pnorm];
    figure(1)
    plot(1:iterations,parray)
    hold on;

    vfeild = vnewfeild;
    ufeild = unewfeild;

      if (pnorm) <tollarance  && n2 > 3
          break
      end 

%    if n2 > 500
%        break
%     end 
    
 end

 %final pressure 
 for i = 1:ygridsize
     for j = 1:xgridsize
            pfinal(i,j) = 0.25*(pfeild(i,j) + pfeild(i,j+1) + pfeild(i+1,j) + pfeild(i+1,j+1));
     end 
 end

 %final u velocity 
  for i = 1:ygridsize
     for j = 1:xgridsize
            ufinal(i,j) = 0.5 * (ufeild(i,j) + ufeild(i+1,j) );
     end 
  end

  %final v velocity 
  for i = 1:ygridsize
     for j = 1:xgridsize
            vfinal(i,j) = 0.5 * (vfeild(i,j) + vfeild(i,j+1) );
     end 
 end
toc
 %% post processing 
 figure(2)
 contourf(pfinal);
 axis equal tight 
 set(gca,'Ydir','reverse')
colormap('jet');
 
  figure(3)
 quiver(ufinal,vfinal);
 set(gca,'Ydir','reverse')
 axis square tight 

 figure(4)
 contourf(ufinal);
 axis equal tight 
 set(gca,'Ydir','reverse')
colormap('jet');




