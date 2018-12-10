clear all
close all
clc

A_nocoup=spconvert(load('A_3d_nocoup.mm'));
A_coup=spconvert(load('A_3d_coup.mm'));
figure
spy(A_nocoup)
figure
spy(A_coup)
C=A_coup-A_nocoup;
figure
spy(C)