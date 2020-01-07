function [posiX, posiY, posiZ] = UeLocation(isd)

% % Set the Seed number
% s = RandStream('mt19937ar','Seed',seed);
% RandStream.setGlobalStream(s);
% % End Setting Seed number

% c = isd/3;
% a = c/2;
% b = a*sqrt(3);
% 
% %posiX = unifrnd(-(c+a), 2*c);
% posiX = unifrnd(-2*(c), 2*c);
% 
% if posiX < -(c*a)
%  posiX = posiX + (c+2*a);
% end
% 
% if posiX >= -c && posiX <= 0
%     posiY = unifrnd(-2*b, 2*b);
% elseif posiX >= a && posiX <= c+a
%     posiY = unifrnd(-b, b);
% elseif posiX < -c
%     bound = (b/a)*(posiX+c+2*a);
%     posiY = unifrnd(2*b-bound, bound);
%     if unifrnd(0,1) >= 0.5
%         posiY = -posiY;
%     end
% elseif posiX > 0 && posiX < a
%     bound = (-b/a)*posiX +2*b;
%     posiY = unifrnd(-bound, bound);
% elseif posiX > c+a
%     bound = (-b/a)*(posiX-2*c);
%     posiY = unifrnd(-bound, bound);
% end

radius = isd/2;
radius = radius * sqrt(rand);
theta = 2 * pi * rand;

posiX = radius * cos(theta);
posiY = radius * sin(theta);

%maxFloor = unidrnd(5,1)+3;
%numFloor = unidrnd(maxFloor,1);
maxFloor = randi(5)+3;
numFloor = randi(maxFloor);
posiZ = 3*(numFloor-1)+1.5;
    
end