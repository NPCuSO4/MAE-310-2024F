function [n_el, n_np, x_coor, y_coor, IEN] = readMsh(filename, n_en)

% Reference: https://gmsh.info/doc/texinfo/#MSH-file-format

fid = fopen(filename, 'r');
line = fgetl(fid); 
while ischar(line)
    if strcmp(line, '$Nodes')
        line = fgetl(fid);
        break
    end
    line = fgetl(fid);
end

num = str2num(line);
n_np = num(2);
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);
IEN = zeros(0, n_en);

cnt = 0;
line = fgetl(fid); 
while ischar(line)
    num = str2num(line);
    if strcmp(line, '$Elements')
        line = fgetl(fid);
        break
    elseif length(num) == 3
        cnt = cnt + 1;
        x_coor(cnt) = num(1);
        y_coor(cnt) = num(2);
    end
    line = fgetl(fid);
end

cnt = 0; 
while ischar(line)
    num = str2num(line);
    if strcmp(line, '$EndElements')
        break
    elseif length(num) == 5
        cnt = cnt + 1;
        tmp = zeros(1,4);
        tmp(1) = num(2);
        tmp(2) = num(3);
        tmp(3) = num(4);
        tmp(4) = num(5);
        IEN = [IEN; tmp];
    end
    line = fgetl(fid);
end
n_el = cnt;

end