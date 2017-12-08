function newind=mutation(pop,params,state,data,i)
%MUTATION    Creates a new individual for GPLAB by mutation.
%   NEWIND=MUTATION(POPULATION,PARAMS,STATE,DATA,PARENT) returns a new
%   individual created by substituting a random subtree of PARENT by a
%   new randomly created tree, with the same depth/size restrictions
%   as the initial random trees.
%
%   Input arguments:
%      POPULATION - the population where the parent is (array)
%      PARAMS - the parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the current dataset for the algorithm to run (array)
%      PARENT - the index of the parent in POPULATION (integer)
%   Output arguments:
%      NEWIND - the newly created individual (struct)
%
%   See also CROSSOVER, APPLYOPERATOR
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

ind=pop(i);

% calculate number of nodes (we need it to pick a random branch)
if isempty(ind.nodes)
   ind.nodes=nodes(ind.tree);
end

x=intrand(1,ind.nodes); % mutation point

% node to mutate (whole branch from this point downwards):
%xnode=findnode(ind.tree,x); 

% build a new branch for this tree, no deeper/bigger than the initial random trees:
% (and obeying the depth/size restrictions imposed by the limits in use)
%newtree=maketree(state.iniclevel,state.functions,state.arity,0,state.depthnodes,xnode.nodeid-1);
newtree=maketree(state.iniclevel,state.functions,state.arity,0,params.depthnodes,[],x-1);
% (the maximum size of the new branch is the same as the initial random trees)
% (0 means no exact level)

% swap old branch with new branch in only one step, as if this were
% crossover (but discard the resulting nind):
nind.tree=newtree;
ind.tree=swapnodes(ind.tree,nind.tree,x,1);
%ind.tree=swapnode(ind.tree,x,newtree);

ind.id=[];
ind.origin='mutation';
ind.parents=[pop(i).id];
ind.xsites=[x];
ind.str=tree2str(ind.tree);
ind.fitness=[];
ind.adjustedfitness=[];
ind.result=[];
ind.testfitness=[];
ind.testadjustedfitness=[];
ind.nodes=ind.tree.nodes;
ind.introns=[];
ind.level=[];

newind=ind;
