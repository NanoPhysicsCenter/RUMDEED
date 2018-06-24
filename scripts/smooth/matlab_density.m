%-------------------------------------------------------------------------%
% Plot density plots for emission and absorbsion                          %
% Kristinn Torfason                                                       %
% 11.05.2018                                                              %
%-------------------------------------------------------------------------%

clear all
close all

% File names for the binary files containing the data 
file_emit = '../../data/Test-Tip/out/density_emit.bin';
file_abs = '../../data/Test-Tip/out/density_absorb_bot.bin';

%--------------------------------------------------------------------------
% Emission
figure()
[x, y, emit] = Read_Density_File(file_emit);

% Filter out data for each emitter
% Emitter number 1
x_1 = x(emit == 1);
y_1 = y(emit == 1);

% Emitter number 2
x_2 = x(emit == 2);
y_2 = y(emit == 2);

% Emitter number 3
x_3 = x(emit == 3);
y_3 = y(emit == 3);

% Emitter number 4
x_4 = x(emit == 4);
y_4 = y(emit == 4);

% Plot the graph
plot(x_1, y_1, 'b.', x_2, y_2, 'r.', x_3, y_3, 'k.', x_4, y_4, 'g.')
title('Emission density scatter')
xlabel('x [nm]')
ylabel('y [nm]')
axis equal
legend('Emitter 1', 'Emitter 2', 'Emitter 3', 'Emitter 4')


figure()
h = histogram2(x, y, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
title('Emission density histogram')
h.NumBins = [100 100];
colorbar
axis equal

%--------------------------------------------------------------------------
% Absorption
figure()
[x, y, emit] = Read_Density_File(file_abs);

% Filter out data for each emitter
% Emitter number 1
x_1 = x(emit == 1);
y_1 = y(emit == 1);

% Emitter number 2
x_2 = x(emit == 2);
y_2 = y(emit == 2);

% Emitter number 3
x_3 = x(emit == 3);
y_3 = y(emit == 3);

% Emitter number 4
x_4 = x(emit == 4);
y_4 = y(emit == 4);

% Plot the graph
plot(x_1, y_1, 'b.', x_2, y_2, 'r.', x_3, y_3, 'k.', x_4, y_4, 'g.')
title('Absorption density scatter')
xlabel('x [nm]')
ylabel('y [nm]')
axis equal
legend('Emitter 1', 'Emitter 2', 'Emitter 3', 'Emitter 4')

figure()
h = histogram2(x, y, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
title('Absorption density histogram')
h.NumBins = [100 100];
colorbar
axis equal

%--------------------------------------------------------------------------
% Read unformated binary stream file from Fortran
% The file has three records for each emission/absorption event
% 1: x coordinate in double precission
% 2: y coordinate in double precission
% 3: number of the emitter wich the particle came from as a 32 bit int
function [x_data, y_data, emit_data] = Read_Density_File(filename)
k = 0;
i = 3;

length_scale = 1.0E-9; % 1 nm

% Open the file for reading
fid = fopen(filename);

%-----------------------------------
% Figure out the length of the file
% store current seek
current_seek = ftell(fid);
% move to end
fseek(fid, 0, 1);
% read end position
file_length = ftell(fid);
% move to previous position
fseek(fid, current_seek, -1);
%----------------------------------

%----------------------------------
% Calculate the number of events recorded in the file.
% We have two double precission values (8 bits) and one 32 bit integer
% (4 bits) for each particle.
N = file_length/(2*8+1*4);
x_data = zeros(1, N);
y_data = zeros(1, N);
emit_data = zeros(1, N);

while (i > 0)
    k = k + 1;
    
    % The first two records are double precision numbers
    [data, i] = fread(fid, 2, 'double');
    if (i == 2) % Check if we actually read some thing
        x_data(k) = data(1) / length_scale;
        y_data(k) = data(2) / length_scale;
    end

    % The third record is a 32bit integer.
    [data, i] = fread(fid, 1, 'int32');
    if (i == 1) % Check if we actually read some thing
        emit_data(k) = data;
    end
end

% Close the file
fclose(fid);
end