function [x, t, y] = lane_change(start, width, slope, type, speed, num, length, varargin)
% Generates the time and coordinates for either a single or double lane change
% manuever at a particular speed.
%
% Parameters
% ----------
% start : float
%   The starting point along the x axis in meters.
% width : float
%   The width of the lane deviation.
% slope : float
%   The slope of the lane change.
% type : string
%   Either 'single' or 'double'. A double lane change return to x = 0.
% speed : float
%   Speed of travel.
% num : integer
%   Number of time steps.
% length : float
%   The length of path.
% laneLength : float, optional
%   Length of the lane for a double lane change.
%
% Returns
% -------
% x : matrix, (num, 1)
%   The longitudinal path.
% y : matrix, (num, 1)
%   The lateral path.
% t : matrix, (num, 1)
%   Time.

x = 0:(length - start) / num:length;
t = x / speed;

y = zeros(length(x), 1);
if strcmp(type, 'single')
    y(start:slope(
