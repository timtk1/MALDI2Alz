%{
This script reads through annotations downloaded from LIPIDMAPS in .tsv
format. Because many possible annotations are generated, this script picks
the closest match in terms of ppm error.
%}


filename = ".\mzlist_pklist3_june22.tsv"; % replace with your file name
data = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');

% Extract unique input masses
unique_masses = unique(data{:, 1});

% Initialize the index array to store the best matches
best_match_indices = zeros(length(unique_masses), 1);

% Loop over each unique mass and find the row with the smallest mass difference
for i = 1:length(unique_masses)
    current_mass = unique_masses(i);
    mass_rows = find(data{:, 1} == current_mass);
    mass_diffs = data{mass_rows, 3}; %calc mass diff
    [~, min_index] = min(mass_diffs); %find indx
    best_match_indices(i) = mass_rows(min_index); %store the best match indices
end

% best indices
best_matches = data(best_match_indices, :);

%write to table or save 
save(best_matches);
%writetable(best_matches, 'best_matches.tsv', 'FileType', 'text', 'Delimiter', '\t');
