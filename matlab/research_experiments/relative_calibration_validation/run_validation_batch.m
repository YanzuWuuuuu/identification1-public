function run_validation_batch(seeds)
% Paired, predeclared runs; no outcome-dependent seed selection.
if nargin<1, seeds=2:8; end
folder=fileparts(mfilename('fullpath')); diary(fullfile(folder,'outputs','batch_log.txt'));
cleanup=onCleanup(@() diary('off'));
for seed=seeds
    for scenario={'common','repeated','distinct'}
        validation_case(scenario{1},seed);
        drawnow;
    end
end
end
