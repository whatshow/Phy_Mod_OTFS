%
% this script tests CE in the EVA model using NMSE
%   > "Evolved universal terrestrial radio access (E-UTRA); base station (BS)
%   radio trans mission and reception, version 8.6.0," document 3GPP TS
%   36.104, Jul. 2009
%
clear;
clc;
%% settings
% OTFS configuration
N = 64;                          % time slot number
M = 16;                          % subcarrier number
% channels
p = 9;
