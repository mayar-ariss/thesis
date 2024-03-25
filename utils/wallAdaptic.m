% Mayar Ariss 25 Mar 2024

function [successOutput, node_name] = wallAdaptic(XPiers_filtered, YPiers_filtered, ZPiers_filtered, XSpandrels_filtered, YSpandrels_filtered, ZSpandrels_filtered, XNodes_filtered, YNodes_filtered, ZNodes_filtered)

% !!!!! look at facade wall only 
% Reshape X, Y, and Z matrices to column vectors
X_nodes = XNodes_filtered(:);
Y_nodes = YNodes_filtered(:);
Z_nodes = ZNodes_filtered(:);

% !!!!! look at facade wall only 
% Reshape X, Y, and Z matrices to column vectors
X_spandrels = XSpandrels_filtered(:);
Y_spandrels = YSpandrels_filtered(:);
Z_spandrels = ZSpandrels_filtered(:);

% !!!!! look at facade wall only 
% Reshape XPiers, YPiers, and ZPiers matrices to column vectors 
X_piers = XPiers_filtered(:, 1:8);
Y_piers = YPiers_filtered(:, 1:8);
Z_piers = ZPiers_filtered(:, 1:8);

X_piers = X_piers(:);
Y_piers = Y_piers(:);
Z_piers = Z_piers(:);

% Determine the number of spandrels and piers
num_nodes= numel(X_nodes);
num_spandrels = numel(X_spandrels);
num_piers = numel(X_piers);
num_elements=num_spandrels+num_piers+num_nodes;

% Create a column vector to indicate whether each element is a spandrel (1) or not (0)
is_spandrel = ones(num_spandrels, 1); % Assuming all elements are spandrels initially

% Create a combined matrix of X, Y, Z, their indices, and the spandrel indicator
spandrels = [X_spandrels, Y_spandrels, Z_spandrels, (1:num_spandrels)', is_spandrel];

% Sort the spandrel combined matrix based on Z first, then Y, and finally on X within each Z-Y group
sorted_spandrels = sortrows(spandrels, [3, 2, 1]);



% Create a column vector to indicate whether each element is a piers (1) or not (0)
is_pier = zeros(num_piers, 1); % Assuming all elements are not piers initially

% Create a combined matrix of XPiers, YPiers, ZPiers, their indices, and the piers indicator
piers = [X_piers, Y_piers, Z_piers, (num_spandrels+1:num_spandrels+num_nodes)', is_pier];

% Sort the piers combined matrix based on Z first, then Y, and finally on X within each Z-Y group
sorted_piers = sortrows(piers, [3, 2, 1]);


% Create a column vector to indicate whether each element is a piers (1) or not (0)
is_node = ones(num_nodes, 1).*2; % Assuming all elements are not piers initially

% Create a combined matrix of XPiers, YPiers, ZPiers, their indices, and the piers indicator
nodes= [X_nodes, Y_nodes, Z_nodes, (num_spandrels+num_nodes+1:num_elements)', is_node];

% Sort the piers combined matrix based on Z first, then Y, and finally on X within each Z-Y group
sorted_nodes = sortrows(nodes, [3, 2, 1]);


% Concatenate the combined matrices for spandrels and piers
combined = [spandrels; piers; sorted_nodes];

% Sort the combined matrix based on Z first, then Y, and finally on X within each Z-Y group
sorted_combined = [sorted_spandrels; sorted_piers; sorted_nodes];

% Extract sorted X, Y, Z coordinates
X_sorted = sorted_combined(:, 1);
Y_sorted = sorted_combined(:, 2);
Z_sorted = sorted_combined(:, 3);

count=2;
Xidx=X_sorted(1);
for i=2:length(X_sorted)
    if X_sorted(i)~=X_sorted(i-1)
        Xidx(count)=X_sorted(i);
        count=count+1;
    end
end

count=2;
Zidx=Z_sorted(1);
for i=2:length(Z_sorted)
    if Z_sorted(i)~=Z_sorted(i-1)
        Zidx(count)=Z_sorted(i);
        count=count+1;
    end
end


% Generate node names based on sorted indices
node_name = cell(size(X_sorted));
x_count=1;
z_count=1;
z_level_prev=1;
x_level_prev=1;

for i = 1:numel(X_sorted)

    % Extracting the indices of z for naming convention
    if Zidx(z_count) == Z_sorted(i)
        z_level = z_level_prev; % Get the z level
    else
        z_count=z_count+1;
        z_level = z_level_prev+1;
        z_level_prev=z_level;
        x_level_prev=0; %indicates change in z value
    end

    % Extracting the indices of x for naming convention
    if Xidx(x_count) == X_sorted(i)
        x_level = x_level_prev; % Get the x level
    else
        x_count=x_count+1;
        x_level = x_level_prev+1;
        x_level_prev=x_level;
    end


    % Generate the node name based on the pattern
    node_name{i} = sprintf('%d%d%d%d%d%d%d', z_level, 0, 0, x_level, 0, 0, i);
end

% Convert node_name to a column vector
node_name = node_name(:);

% Convert numeric arrays to cell arrays
X_sorted_cell = num2cell(X_sorted);
Y_sorted_cell = num2cell(Y_sorted);
Z_sorted_cell = num2cell(Z_sorted);
index_cell= num2cell(sorted_combined(:,4));
s_or_p_cell = num2cell(sorted_combined(:, 5));

% Combine sorted data with node names
sorted_data_with_names = [X_sorted_cell, Y_sorted_cell, Z_sorted_cell, index_cell, s_or_p_cell, node_name];

% Match node names with original combined matrix indices
combined_with_names = cell(size(combined, 1), size(combined, 2) + 1); % Initialize combined_with_names
combined_with_names(:, 1:end-1) = num2cell(combined); % Copy X, Y, Z, indices to combined_with_names

% Match node names with corresponding indices
for i = 1:size(combined, 1)
    idx = find(sorted_combined(:, 4) == i, 1); % Get the index in the sorted_combined matrix
    combined_with_names{i, end} = node_name{idx}; % Assign the corresponding node name
end

% Define the number of nodes per element
nodes_per_element = 8;

% Determine the number of elements
num_elements = size(combined_with_names, 1) / nodes_per_element;

% Initialize the matrix to store elements and their node names
elements_with_nodes = cell(num_elements, 1);

% Extract node names for each element
for i = 1:num_elements
    % Determine the indices for the current element
    start_idx = (i - 1) * nodes_per_element + 1;
    end_idx = start_idx + nodes_per_element - 1;

    % Extract node names for the current element
    node_names = combined_with_names(start_idx:end_idx, end);

    % Store node names for the current element
    elements_with_nodes{i} = node_names';
end

% Display the elements with their node names
for i = 1:num_elements
    fprintf('Element %d: %s\n', i, strjoin(elements_with_nodes{i}, ' '));
end


% Save X_sorted, Y_sorted, Z_sorted, and node_name to a text file on the desktop
output_file = '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/sorted_coordinates.txt';
fid = fopen(output_file, 'w');
if fid == -1
    error('Could not open file for writing: %s', output_file);
end

fprintf(fid, 'nod.name\tx\ty\tz\n');
for i = 1:numel(X_sorted)
    fprintf(fid, '%s\t%f\t%f\t%f\n', node_name{i}, X_sorted(i), Y_sorted(i), Z_sorted(i));
end


% Save elements_with_nodes to a text file on the desktop
output_file = '/Users/mayar/Desktop/ICL/4Y/Thesis/Eddy/elem_coordinates.txt';
fid = fopen(output_file, 'w');
if fid == -1
    error('Could not open file for writing: %s', output_file);
end

fprintf(fid, 'elm.name\n');
for i = 1:numel(elements_with_nodes)
    % Concatenate node names of the current element with spaces
    node_names_str = strjoin(elements_with_nodes{i}, ' ');

    % Print the concatenated node names
    fprintf(fid, '%s\n', node_names_str);
end


fclose(fid);
successOutput = true;

end
