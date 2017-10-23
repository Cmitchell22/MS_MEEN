clc; clear all; close all;
syms 't';

alpha = 0.5;
fdefun = t^2;
t0 = 1;
tfinal = 100;
y0 = 0;
h = 1;
fde_BongsooJang(alpha,fdefun,t0,tfinal,y0,h);