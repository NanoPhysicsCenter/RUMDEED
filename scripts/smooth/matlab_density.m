% Read unformated stream file from Fortran

clear all
close all

file_emit = '../../data/FE-Test/out/density_emit.bin';
file_abs = '../../data/FE-Test/out/density_absorb.bin';

%--------------------------------------------------------------------------
% Emission
figure(1)
[x, y, emit] = Read_Density_File(file_emit);

plot(x, y, 'b.')
title('Emission density')
xlabel('x [nm]')
ylabel('y [nm]')

%--------------------------------------------------------------------------
% Absorption
figure(2)
[x, y, emit] = Read_Density_File(file_abs);

plot(x, y, 'b.')
title('Absorption density')
xlabel('x [nm]')
ylabel('y [nm]')

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

fclose(fid);
end
