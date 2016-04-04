function index = getStructIndex(origStruct, searchStruct)
% Returns indices of searched substructure in cell array of structures.
%
% Input:
%   origStruct   - cell array of structures for search
%   searchStruct - substructure to be searched in 'origStruct'
%
% Output:
%   index - vector of indices of searched substructure in 'origStruct'

  index = [];
  searchedFields = findFields(searchStruct);
  nFields = length(searchedFields);
  searchedValues = cell(1, nFields);
  for j = 1:nFields % find all field values
    searchedValues{j} = eval(['searchStruct.', searchedFields{j}]);
  end
  for i = 1:length(origStruct)
    correctFields = true;
    for j = 1:nFields % compare all needed fields
      correctFields = correctFields && all(isequal(eval(['origStruct{i}.', searchedFields{j}]), searchedValues{j}));
    end
    if correctFields
      index(end+1) = i;
    end
  end
end

function resFields = findFields(str)
% find recursive fields in structure

  resFields = {};
  strFields = fieldnames(str);
  nFields = length(strFields);
  for i = 1:nFields
    if isstruct(str.(strFields{i}))
      actualFields = findFields(str.(strFields{i}));
      resFields = [resFields, cellfun(@(x) [strFields{i}, '.', x], actualFields, 'UniformOutput', false)];
    else
      resFields = [resFields, strFields{i}];
    end
  end
end