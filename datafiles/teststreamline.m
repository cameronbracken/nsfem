function teststreamline

load wind
[sx,sy,sz] = meshgrid(80,20:10:50,0:5:15);

% Plot the streamlines
figure(1)
h = streamline(x,y,z,u,v,w,sx,sy,sz);
set(h,'Color','red')
view(3)
title('Streamline Plot - See how fast I am')

% Plot the streamlines with color corresponding to magnitude of vector
% field
figure(2)
c = sqrt(u.^2 + v.^2 + w.^2);
streamlinec(stream3(x, y, z, u, v, w, sx, sy, sz), x, y, z, c)
title('Colored Streamline Plot - Boy Am I Slow!')
view(3)

return

function streamlinec(s, X, Y, Z, C)

% Iterate over streamlines
for slines = 1:length(s)

    current_streamline = s{slines};

    % Interpolate the color field onto the current streamline
    ci = interp3(X, Y, Z, C, current_streamline(:,1), ...
        current_streamline(:,2), current_streamline(:,3));

    % Plot each segment of the streamline as a patch object, color the
    % segment according to the ci vector
    for seg = 1:length(current_streamline)-1
        patch('vertices', current_streamline, 'faces', seg:seg+1, ...
            'facevertexcdata', ci, ...
            'Facecolor', 'none', 'edgecolor', 'flat')
    end
end

return 