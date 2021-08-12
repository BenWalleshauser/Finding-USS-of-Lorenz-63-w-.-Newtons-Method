%Ben Walleshauser
%7/19

%Finding the USS of Lorenz '63 with Newton's Method
%Could create basins of where initial conditions land 
%Some cool numerical stuff: https://www.youtube.com/c/ProfJeffreyChasnov/videos
%% Initializing
clear
clc
close all

FirstGuess = [-20; -20; 20];
NumIts = 10;    %Doesn't require 'a lot' of iterations to get a really good estimate

sigma = 10;
rho = 28;
beta = 8/3;

%% Newton's Method 

u = zeros(3, NumIts+1);
u(:,1) = FirstGuess;

for i = 1:NumIts
    x = u(1,i);
    y = u(2,i);
    z = u(3,i);
    u(:,i+1) = u(:,i) - ([-rho rho 0; rho-z -1 -x; y x -beta]^-1)*[sigma*(y-x); x*(rho-z)-y; x*y-beta*z];
end

format long
USS_predicted = u(:,NumIts+1)
%% Checking

%Analytic values:
USS1 = [0; 0; 0]                                                        
USS2 = [sqrt(beta*(rho-1)); sqrt(beta*(rho-1)); rho-1]
USS3 = [-sqrt(beta*(rho-1)); -sqrt(beta*(rho-1)); rho-1]

%% Nullclines

[X1,Z1] = meshgrid(-30:30,-120:4:120);
Y1 = X1;
%surf(X1,Y1,Z1)

[X2,Y2] = meshgrid(-30:1:30,-30:1:30);
Z2 = rho-Y2./X2;
%surf(X2,Y2,Z2)

[X3,Y3] = meshgrid(-30:30,-30:30);
Z3 = X3.*Y3./beta;
%surf(X3,Y3,Z3)

[t, X] = ode45('lorenz', [0:0.01:35], [10 10 10]); 
hold on
s.EdgeColor = 'none';
surf(X1,Y1,Z1,'facecolor','b','FaceAlpha',0.5,'EdgeColor','none')
surf(X2,Y2,Z2,'facecolor','r','FaceAlpha',0.5,'EdgeColor','none')
surf(X3,Y3,Z3,'facecolor','g','FaceAlpha',0.5,'EdgeColor','none')
plot3(X(:,1),X(:,2),X(:,3))
plot3(USS1(1),USS1(2),USS1(3),'o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
plot3(USS2(1),USS2(2),USS2(3),'o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
plot3(USS3(1),USS3(2),USS3(3),'o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
legend('First Nullcline','Second Nullcline','Third Nullcline','Lorenz63','USS1','USS2','USS3')
xlabel('x')
ylabel('y')
zlabel('z')
title('USS of Lorenz 63')
hold off


%% Basins of USS

%Won't be easy to visualize in R3 but could do a cross section 'movie'











