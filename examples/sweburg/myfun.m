function F = myfun(v)
global x t
F = v-sin(2*pi*(x-v*t))*0.4;