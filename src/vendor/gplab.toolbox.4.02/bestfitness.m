function [validind,state]=bestfitness(ind,pop,params,state,data,parentindices)
%STRICTDEPTH    Applies strict depth filters to a new GPLAB individual.
%   [VALIDIND,STATE]=STRICTDEPTH(IND,POP,PARAMS,STATE,DATA,PARENTS)
%   tests if an individual (IND) conforms to the strict maximum
%   depth rules. If not, returns one of its parents.
%
%   Input arguments:
%      IND - individual to be validated (array)
%      POPULATION - the current population of the algorithm (array)
%      PARAMS - the running parameters of the algorithm (struct)
%      STATE - the current state of the algorithm (struct)
%      DATA - the dataset for use in the algorithm (struct)
%      PARENTS - the indices of the parents of IND (matrix)
%   Output arguments:
%      VALIDIND - valid individual (IND or one of its parents) (array)
%      STATE - the updated state of the algorithm (struct)
%
%   References:
%      Koza J. (1992) "Genetic Programming: On the Programming of
%      Computers by Means of Natural Selection". MIT Press.
%      Silva S. and Almeida J. (2003) "Dynamic maximum tree depth - a
%      simple technique for avoiding bloat in tree-based GP". GECCO-2003.
%
%   See also VALIDATEINDS, STRICTNODES, HEAVYDYNDEPTH, ... (the other filters)
%
%   Copyright (C) 2003-2015 Sara Silva (sara@fc.ul.pt)
%   This file is part of the GPLAB Toolbox

% we need fitness:
[ind,state]=calcfitness(ind,params,data,state,0);

if or(and(params.lowerisbetter,ind.fitness<=pop(parentindices(1)).fitness),and(~params.lowerisbetter,ind.fitness>=pop(parentindices(1)).fitness))
    state.lastid=state.lastid+1;
    ind.id=state.lastid;
else
      % else, substitute by its parent:
      ind=pop(parentindices(1));
      % tree nodes is needed:
      if isempty(ind.nodes)
		ind.nodes=treenodes(ind.tree);
      end
      if length(state.rejected)<ind.nodes
	state.rejected(ind.nodes)=1;
      else
	state.rejected(ind.nodes)=state.rejected(ind.nodes)+1;
      end
end       

validind=ind;
