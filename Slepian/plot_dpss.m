clc;
clear;
close all;

fileID = fopen('DPSS.bin','rb');
A = fread(fileID,[64,4],'double');
plot(A);
xlabel('n');
title('Discrete Prolate Spheroidal Sequences');
%hold on;
seq = dpss(64,4,4);
%plot(seq);
figure;
plot(seq-A);
title('Difference');

