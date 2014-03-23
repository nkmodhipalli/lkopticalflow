% Copyright 2014 Naresh Kumar Modhipalli
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function V = getlkflow(I,J,U,lkwin,levels)

k = [1,4,6,4,1;4,16,24,16,4;6,24,36,24,6;4,16,24,16,4;1,4,6,4,1;]/256;
kd = [0.5 0 -0.5];

lkwin2 = floor(lkwin/2);
num_pts = size(U,2);
V = zeros(size(U));

pyrI = cell(levels,1);
pyrJ = cell(levels,1);
pyrI{1} = I;
pyrJ{1} = J;
for l = 2:levels
    pyrI{l} = conv2(pyrI{l-1},k,'same');
    pyrI{l} = pyrI{l}(1:2:end,1:2:end);
    pyrJ{l} = conv2(pyrJ{l-1},k,'same');
    pyrJ{l} = pyrJ{l}(1:2:end,1:2:end);
end

pyrIGradX = cell(levels,1);
pyrIGradY = cell(levels,1);
for l = 1:levels
    pyrIGradX{l} = conv2(pyrI{l},kd,'same');
    pyrIGradY{l} = conv2(pyrI{l},kd','same');
end

for pt = 1:num_pts
    u = U(:,pt);
    g = [0 0]';
    for l = levels:-1:1
        u_l = ceil(u/2^(l-1));
        G(1,1) = sum(sum((pyrIGradX{l}(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
            max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2))).^2));
        G(2,2) = sum(sum((pyrIGradY{l}(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
             max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2))).^2));
        G(2,1) = sum(sum((pyrIGradX{l}(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
            max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2))).*...
            (pyrIGradY{l}(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
            max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2)))));
        G(1,2) = G(2,1);
        v = [0 0]';
        for iter = 1:20
            del_ind = v+g;
            tempJ = zeros(size(pyrJ{l}));
            tempJ(1:end - abs(del_ind(1)),1:end - abs(del_ind(2))) = pyrJ{l}(max(1,1+del_ind(1)):min(end,end+del_ind(1)),...
                max(1,1+del_ind(2)):min(end,end+del_ind(2)));
            delI = (pyrI{l}(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
                max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2))) - ...
                (tempJ(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
                max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2)));
            b(1,1) = sum(sum(delI.*(pyrIGradX{l}(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
                max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2)))));
            b(2,1) = sum(sum(delI.*(pyrIGradY{l}(max(1,u_l(1)-lkwin2):min(end,u_l(1)+lkwin2),...
                max(1,u_l(2)-lkwin2):min(end,u_l(2)+lkwin2)))));
            eta = flipud(pinv(G)*b);
            v = round(v + eta);
        end
        d = v;
        temp = g+d;
        g = 2*(g+d);
    end
    d = temp;
    v = u+d;
    V(:,pt) = v;
end
