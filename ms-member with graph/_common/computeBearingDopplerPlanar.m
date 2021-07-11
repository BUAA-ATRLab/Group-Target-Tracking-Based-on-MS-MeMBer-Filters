function [bearing, frequency] = computeBearingDopplerPlanar(sensorPos, targetState, lambda)
%   Measurement function to compute the bearing and Doppler-shift frequency
%   given sensor and target positions and target frequency
%   We can have multi-sensor to single-target or vice versa, but we can not 
%   have multi-sensor to multi-target
%
%   Input: 
%       sensorPos: 2xN matrix where column row is one sensor's position where
%       N is number of sensors 
%       targetState: 4xJ matrix where each column is one target's state where
%       J is number of targets (J=1 if N~=1. If J~=1, N=1)
%       The target state is [x, dot(x), y, dot(y)]', where dot() indicates
%       speed
%       lambda: 1xN vector, each entry contains the signal wavelength for
%       the corresponding frequency measurement
%
%   Output:
%       bearing: single row vector of target bearing relative to sensor in
%       degrees
%       frequency: single row vector of Doppler-shift frequency in Hz
%       Both bearing and fequency are of the format 
% snsor     ____1_______2_______._______._______N____
% targets   | 1...J | 1...J | 1...J | 1...J | 1...J |      
%          
%% Input validation
N = size(sensorPos,2); % number of sensors
J = size(targetState,2); % number of targets

% First, we have to make sure we have the right number of sensors/targets
if (N~=1&&J~=1)
    error('One sensor and multiple targets or multiple sensors and one target');
end

targetPos = targetState([1,2],:); 
targetVelocity = targetState([3,4],:);

% Replicate sensor/target positions so they have the same length
if (N==1)
    sensorPos = repmat(sensorPos, 1, J);
else
    if (J==1)
        targetPos = repmat(targetPos, 1, N);
        targetVelocity = repmat(targetVelocity, 1, N);
    end
end

%% Compute measurements
% center coordinate system in the sensor
difX = targetPos(1,:) - sensorPos(1,:);
difY = targetPos(2,:) - sensorPos(2,:); 
% Compute the bearing in degrees
bearing = atan2(difY, difX)/pi*180;  % output bearing in degrees 

% Compute the Doppler-shift frequency
frequency = zeros(size(bearing));
if sum(difX+difY) ~= 0
    frequency = 2./lambda.*((difX.*targetVelocity(1,:)+difY.*targetVelocity(2,:))./sqrt(difX.^2+difY.^2));
end

end
