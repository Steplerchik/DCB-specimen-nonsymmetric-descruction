clc
clear all
close all
syms t r

U=2+r*t;
r=solve(U==0,r)
%dU=diff(U, r)
%d2U=diff(dU, t)

t=4;
r=eval(r)