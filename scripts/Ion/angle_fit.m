clear all
close all

energy = importdata('data_e.txt');
angle = importdata('data_angle.txt');

plot(energy, angle)