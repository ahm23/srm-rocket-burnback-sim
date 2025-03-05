function [field] = genDistanceField(r_outer, geometry)
% GENDISTANCEFIELD Generates a minimum distance field for the fuel grain.
%   Resolution options can be adjusted in `configure.m`
% Inputs:
%    r_outer    : outer radius of the fuel grain
% Outputs:
%    field      : a field of points

% Define initial grain geometry manually (tube cross-section)
theta = linspace(0, 2*pi, r_outer_resolution);
x_outer = r_outer * cos(theta);
y_outer = r_outer * sin(theta);
x_inner = geometry(:,1);
y_inner = geometry(:,2);

polygon = [x_outer, x_inner; y_outer, y_inner]';

%% Create Initial Field Grid
[X, Y] = meshgrid(min(x_outer)-2:grid_size:max(x_outer)+2, min(y_outer)-2:grid_size:max(y_outer)+2);
points = [X(:), Y(:)];

%% Minimum Distance Function
field = inf(size(X));
for i = 1:length(x_inner)
    dist = sqrt((X - x_inner(i)).^2 + (Y - y_inner(i)).^2);
    field = min(field, dist);
end

%% Point-In-Polygon Algorithm
[in, on] = inpoly2(points, polygon);
positive_points = points(in & ~on, :);  % Points on the fuel grain
boundary_points = points(on, :);        % Points on the boundary (outer/inner)
outside_points = points(~in, :);
[in, on] = inpoly2(outside_points, polygon_inner);
negative_points = outside_points(in & ~on, :);  % Points in the grain cavity
outside_points = outside_points(~in, :);        % Points outside the casing

% Make `negative_points` negative in `field`
for i = 1:size(negative_points,1)
    % Find the corresponding indices in the grid
    idx = (X == negative_points(i,1)) & (Y == negative_points(i,2));
    field(idx) = -field(idx); % Multiply by -1
end

% Make `outside_points` NaN in `field`
for i = 1:size(outside_points,1)
    % Find the corresponding indices in the grid
    idx = (X == outside_points(i,1)) & (Y == outside_points(i,2));
    field(idx) = NaN;
end

end