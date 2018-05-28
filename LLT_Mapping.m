function [Av,Bv,target]=GLT_Mapping(Lv, Ln)
% This program recovers the blurred channel via the sharp channel of exact
% same scene (example: same depth).
% 'Lv' denotes the blurred one while the 'Ln' denotes the sharp one.
% 'Av' & 'Bv' denotes parameters matices.

[m, n, ~] = size(Lv);

blurredLv = imgaussfilt(Lv, 10);
blurredLn = imgaussfilt(Ln, 10);

% Initialization
Av =  0.8*ones(m, n);
Bv =  0.05*ones(m, n);
[Av Bv] = gradDes(Av, Bv, blurredLv, blurredLn);

target=Av.*stackSharp+Bv;

