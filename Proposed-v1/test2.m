clc;clear;
a = randn(1024,1) + 1j * randn(1024,1);
tic
    fft(a,1024);
toc

delta = [1:1024].';
F = exp(1j * 2 * pi * delta *delta.' / 1024);

tic
    F * a;
toc