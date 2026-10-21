function [vecM,C2,B2]=get_truth()
% Returns rigid-body vec(M) and flexible mode-shape matrices from 3/model.m
here=fileparts(mfilename('fullpath')); root=fileparts(here);
old=pwd; c=onCleanup(@() cd(old)); cd(fullfile(root,'3')); no=4; ni=4; model; %#ok<NASGU>
end
