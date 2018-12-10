close all
clear all
clc

A=load('A_fede.mm');
B=load('Btt_fede.mm');
S=spconvert(A);
C=spconvert(B);
figure
spy(S)
figure
spy(C)