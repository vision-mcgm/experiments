function [ output_args ] = genParams( params )
%Generates block params for experiment, and runs

%Dummy params

%Here we vary target length across blocks, sample length within.

params(1).name='tl';
params(1).type='enum';
params(1).list=[10 25 50];
params(1).scheme='acrossblock';


params(2).name='ratio';
params(2).type='enum';
params(2).list=[1.2 1.4 1.6 1.8 2];
params(2).scheme='withinblock';


params(3).name='wait';
params(3).type='enum';
params(3).list=[1 4 8 10];
params(3).scheme='withinblock';

blockReps=2; %Number of repetitions of each style of block
conditionReps=100; %Number of repetitions of each condition. Must be divisible by blockReps
conditionRepsPerBlock=conditionReps/blockReps;

if mod(blockReps,conditionReps)
    error('conditionReps is not divisble by blockReps');
end


%Get cardinalities
nFactors=size(params,2);
nConditions=1;
for i=1:nFactors
    params(i).n=size(params(i).list,2);
    nConditions=nConditions*params(i).n;
end


%Work out block structure
nBlocks=nConditions*blockReps;
nTrialsPerBlock=

trials=zeros(nBlocks,nTrialsPerBlock,nFactors);

%Work out block params
for f=1:nFactors
    if strcmp(params(f),'acrossblock')
        block

for b=1:nBlocks
    
    
%Run



end

