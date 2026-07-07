# IMPORTS

import numpy as np
import pandas as pd

import os

import math as m
import statistics as stat

import matplotlib
# matplotlib.use('Qt5Agg') 
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import seaborn as sns
sns.set_theme()

# Memory map file
script_dir = os.path.dirname(os.path.abspath(__file__))
lp_filename = os.path.join(script_dir, 'out', 'lp_field_0.bin')
lp_dt_abs_type = np.dtype([('i', np.int64), ('j', np.int64), ('field', np.float64), ('emit',np.int64)])
lp_data_mem_abs = np.memmap(lp_filename, dtype=lp_dt_abs_type, mode='r', order='F')
lp_data = pd.DataFrame.from_records(data=lp_data_mem_abs, columns=lp_data_mem_abs.dtype.names)

script_dir = os.path.dirname(os.path.abspath(__file__))
ic_filename = os.path.join(script_dir, 'out', 'ic_field_0.bin')
ic_dt_abs_type = np.dtype([('i', np.int64), ('j', np.int64), ('field', np.float64), ('emit',np.int64)])
ic_data_mem_abs = np.memmap(ic_filename, dtype=ic_dt_abs_type, mode='r', order='F')
ic_data = pd.DataFrame.from_records(data=ic_data_mem_abs, columns=ic_data_mem_abs.dtype.names)

rms = 0
relative_rms = 0
Nx = 200; Ny = 200
Lx = 10; Ly = 10
lp_field = np.zeros((Nx,Ny))
ic_field = np.zeros((Nx,Ny))
diff_field = np.zeros((Nx,Ny))
x = np.linspace(-Lx,Lx,Nx)
y = np.linspace(-Lx,Ly,Ny)
X, Y = np.meshgrid(x,y)
X = X.T; Y = Y.T
for k in range(Nx*Ny):
    i = lp_data['i'][k]
    j = lp_data['j'][k]

    lp_field[i][j] = lp_data['field'][k]*1e-6
    ic_field[i][j] = ic_data['field'][k]*1e-6
    diff_field[i][j] = (ic_field[i][j] - lp_field[i][j])/abs(lp_field[i][j])*100

    rms += abs(lp_field[i][j]-ic_field[i][j])
    # print(i,j,lp_field[i][j])

rms = rms/(Nx*Ny)

print(f'Laplace average field: {lp_field.mean():.3f} MV/m')
print(f'Image charge average field: {ic_field.mean():.3f} MV/m')
rmse = np.linalg.norm(lp_field-ic_field)/np.sqrt(Nx*Ny)
print(f'RMSE:{rmse:.3f} MV/m')
relative_rms = np.linalg.norm(lp_field-ic_field)/np.linalg.norm(lp_field)
print(f'Relative RMSE: {relative_rms:.3f} MV/m')

rotate = -75
elevation = 20

fig = plt.figure(figsize=(16,8))
ax1 = fig.add_subplot(111, projection='3d', facecolor='white')  # Set facecolor to white to remove grey background
ax1.view_init(elev=elevation, azim=rotate)  # Adjust elevation and azimuth for rotation
ax1.plot_surface(X,Y,lp_field,vmin=lp_field.min(),vmax=lp_field.max(),cmap='viridis',linewidth=0.5, antialiased=True)
# ax1.set_title('Laplace solution')
ax1.set_xlabel('x [nm]',size='xx-large', labelpad=10)
ax1.set_ylabel('y [nm]',size='xx-large', labelpad=10)
ax1.set_zlabel('Field [mV/nm]',size='xx-large', labelpad=20)
ax1.tick_params(axis='both', labelsize='xx-large')
plt.show()
fig = plt.figure(figsize=(16,8))
ax2 = fig.add_subplot(111, projection='3d', facecolor='white') 
ax2.view_init(elev=elevation, azim=rotate)
ax2.plot_surface(X,Y,ic_field,vmin=ic_field.min(),vmax=ic_field.max(),cmap='viridis',linewidth=0.5, antialiased=True)
# ax2.set_title('Image charge solution')
ax2.set_xlabel('x [nm]',size='xx-large', labelpad=10)
ax2.set_ylabel('y [nm]',size='xx-large', labelpad=10)
ax2.set_zlabel('Field [mV/nm]',size='xx-large', labelpad=20)
ax2.tick_params(axis='both', labelsize='xx-large')
plt.show()
fig = plt.figure(figsize=(16,8))
ax3 = fig.add_subplot(111, projection='3d', facecolor='white')
ax3.view_init(elev=elevation, azim=rotate)
ax3.plot_surface(X,Y,diff_field,vmin=diff_field.min(),vmax=diff_field.max(),cmap='viridis',linewidth=0.5, antialiased=True)
# ax3.set_title('Difference between Laplace and image charge solution')
ax3.set_xlabel('x [nm]',size='xx-large', labelpad=10)
ax3.set_ylabel('y [nm]',size='xx-large', labelpad=10)
ax3.set_zlabel('Relative difference [%]',size='xx-large', labelpad=20)
ax3.tick_params(axis='both', labelsize='xx-large')
plt.show()

# x_index = Nx//2
# y_index = Ny//2

# lp_x_data = [lp_field[x_index][j] for j in range(Ny)]
# ic_x_data = [ic_field[x_index][j] for j in range(Ny)]
# lp_y_data = [lp_field[i][y_index] for i in range(Nx)]
# ic_y_data = [ic_field[i][y_index] for i in range(Nx)]

# fig = plt.figure(figsize=(16,8))
# ax1 = fig.add_subplot(211)
# ax1.set_title(f'Laplace vs image charge solution along x = {x_index}')
# ax1.plot(x,lp_x_data,label='Laplace')
# ax1.plot(x,ic_x_data,label='Image charge')
# ax1.set_xlabel('x')
# ax1.set_ylabel('Field [MV/m]')
# ax1.legend()
# ax2 = fig.add_subplot(212)
# ax2.set_title(f'Laplace vs image charge solution along y = {y_index}')
# ax2.plot(y,lp_y_data,label='Laplace')
# ax2.plot(y,ic_y_data,label='Image charge')
# ax2.set_xlabel('x')
# ax2.set_ylabel('Field [MV/m]')
# ax2.legend()

# plt.tight_layout()
# plt.show()
