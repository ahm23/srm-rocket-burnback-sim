function [points] = finocyl_pointy(radius, fins, fin_length, fin_span)
%FINOCYL_POINTY Generate 1D finocyl profile with pointed fins
%   TODO: docstring

%% Precison Parameters
N_arc = 20;         % Points on each arc
N_fin_seg = 10;     % Points on each straight fin segment
tol = 1e-12;        % Tolerance for duplicate point checks

% Fin center angles
theta_c = linspace(0, 2*pi, fins+1);    % NOTE: +1 to generate N sections
theta_c(end) = []                       % NOTE: last is 2pi which overlaps 0, so it is removed

% Initialize an empty array to store the (x,y) points
points = [];

% Helper function for wrapping angles into [0, 2*pi)
wrap2pi = @(a) mod(a, 2*pi);

%% Loop over each fin and generate segments of the profile
prev_fin_end = theta_c(end) + fin_span/2;  % Start from last fin's end

for i = 1:fins
    current_center = theta_c(i);
    fin_start = current_center - fin_span/2;
    fin_end   = current_center + fin_span/2;
    
    % Circle arc from previous fin's end to current fin's start
    a1 = wrap2pi(prev_fin_end);
    a2 = wrap2pi(fin_start);
    if a2 < a1, a2 = a2 + 2*pi; end
    arc_angles = linspace(a1, a2, N_arc)';
    arc_seg = [radius*cos(arc_angles), radius*sin(arc_angles)];
    
    for j = 1:size(arc_seg,1)
        new_pt = arc_seg(j,:);
        if isempty(points) || ~ismembertol(new_pt, points, tol, 'ByRows', true)
            points = [points; new_pt];
        end
    end
    
    % Fin segment from circle at fin start to fin tip
    x_fs = radius*cos(fin_start);
    y_fs = radius*sin(fin_start);
    x_tip = (radius + fin_length)*cos(current_center);
    y_tip = (radius + fin_length)*sin(current_center);
    t = linspace(0,1,N_fin_seg)';
    seg1 = [(1-t)*x_fs + t*x_tip, (1-t)*y_fs + t*y_tip];
    
    for j = 1:size(seg1,1)
        new_pt = seg1(j,:);
        if isempty(points) || ~ismembertol(new_pt, points, tol, 'ByRows', true)
            points = [points; new_pt];
        end
    end
    
    % Fin segment from fin tip back to circle at fin end
    x_fe = radius*cos(fin_end);
    y_fe = radius*sin(fin_end);
    seg2 = [(1-t)*x_tip + t*x_fe, (1-t)*y_tip + t*y_fe];
    
    for j = 1:size(seg2,1)
        new_pt = seg2(j,:);
        if isempty(points) || ~ismembertol(new_pt, points, tol, 'ByRows', true)
            points = [points; new_pt];
        end
    end
    
    % Prepare for next iteration
    prev_fin_end = fin_end;
end

%% Final circle arc from the last fin's end to the first fin's start
first_fin_start = theta_c(1) - fin_span/2;
a1 = wrap2pi(prev_fin_end);
a2 = wrap2pi(first_fin_start);
if a2 < a1, a2 = a2 + 2*pi; end
arc_angles = linspace(a1, a2, N_arc)';
arc_seg = [radius*cos(arc_angles), radius*sin(arc_angles)];

for j = 1:size(arc_seg,1)
    new_pt = arc_seg(j,:);
    if isempty(points) || ~ismembertol(new_pt, points, tol, 'ByRows', true)
        points = [points; new_pt];
    end
end

%% Forced Line Closure Check
% If final point and first point are not the same, append the first point
if norm(points(1,:) - points(end,:)) > tol
    points = [points; points(1,:)];
end
end

