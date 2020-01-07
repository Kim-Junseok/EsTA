function R=GaussRWconfinedCircleReflectingBoundaries(m,d,sigma,X0,Y0)
% Generates a random walk (RW) with m steps.
% The x and y components of each step are drawn from a normal distribution with standard deviation sigma.
% The RW is confined to a circle of diameter d with reflecting boundary conditions
% X0 and Y0 are the initial coordinates of the random walk inside the
% circle. Output is a vector of positions of the particle.
%
%   Usage Examples:
%   GaussRWconfinedCircleReflectingBoundaries(100,50,5,0,0);
%
%   Routine written by Sanaz Sadegh, Colorado State University.
%   Version 1.00
%   December, 2016
xn=randn(m,1); %standard normal r.v.
yn=randn(m,1); %standard normal r.v.
%xn=ones(m,1);
%yn=ones(m,1);
Dx=sigma*xn; %Gaussian(0,sigma) r.v.
Dy=sigma*yn; %Gaussian(0,sigma) r.v.
M=[Dx Dy];
R=zeros(m+1,2);
old=[X0 Y0]; % start RW at X0 and Y0
R(1,:)=old;
for k=1:m
    new=old+M(k,:); % new point
    
    while norm(new)>d/2 % If new point is outside circle use reflective boundar condition untill new point is inside
        coefficients = polyfit([new(1), old(1)], [new(2), old(2)], 1);  % line passing old point and new point which is outside circle
        slope = coefficients (1);
        intercpt = coefficients (2);
        [xout,yout] = linecirc(slope,intercpt,0,0,d/2); % find intrsection of the line with circle
        d1=sqrt((xout(1)-new(1))^2+(yout(1)-new(2))^2); %There are two intersections so calculate the distance from
        d2=sqrt((xout(2)-new(1))^2+(yout(2)-new(2))^2); %new point to these intersection and choose the closer one
        if d1<d2
            P(1)=xout(1);
            P(2)=yout(1);
        else
            P(1)=xout(2);
            P(2)=yout(2);
        end
        a=new-P; %Find the vector from the intersection point to the new point 
        new=new-2*(dot(P,a)/norm(P)^2)*P; % Find the reflection of vector a with respect to the tangent line at the intersection point
        old=P;
    end
    
    old=new;
    R(k+1,:)=old;
    
end
figure(1);
axis([-d/2 d/2 -d/2 d/2])
hold on
plot(X0,Y0,'Marker','o')
plot(R(:,1), R(:,2));
circle([0,0],d/2,1000,'--');
hold off
