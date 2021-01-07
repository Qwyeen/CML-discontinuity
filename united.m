% Unit of time: ms

tProb = 2000; % theorectical max interval between 2 firing event
p0 = 0.1; % firing probability with no influence of previous firing
t0 = 400; % absolute refractory period
t1 = 500; % absolute + relative refractory period

[prob, TimeProb] = ProbabilityFunction(t0, t1, tProb, p0);

figure (1)
plot(TimeProb, prob);
ylim([0,1]);
xlabel('time(ms)');
ylabel('probability');
title('Probability Function (linear increase)');

FN_x = 250; % width of a FN segment
FN_y = 1000; % length of a FN segment
gap = 50; % length of gap
num_segment = 2; % number of FN segments
FN_rgb = [200,200,200]/256; % FN color
[FN_x_auto, FN_y_auto] = FN_map(FN_x, FN_y, gap, num_segment);
figure (2)
patch(FN_x_auto,FN_y_auto,FN_rgb); % place the FN map
axis equal
title('Fibronectin Map');

section_w = 25; % width of FN section
cell_max_x = 120; % maximum cell x direction
cell_min_x = 50; % minimum cell x direction
[cell_info,cell_num, cell_center] = cells(FN_x, FN_y, gap, num_segment, section_w, cell_max_x, cell_min_x);
figure (3)
title('Fibronectin Map with Confluent Cells');
patch(FN_x_auto,FN_y_auto,FN_rgb); % place the FN map
axis equal
hold on
for i = 1 : cell_num
    rectangle('Position',cell_info{i}(1:4))
    hold on
end

[distance_matrix] = dist_mat_calculator(cell_center, cell_num);
sim_time = 2000; % simulation time
cells_lt = zeros(cell_num, sim_time); % local time matrix for all cells
flag = zeros(cell_num,1); % records the time of last firing for cells
thres_dict = 50;

for tt = 1 : sim_time % loop through time
    firingEvent = rand(cell_num,1); % create a random number vector
    P = prob(cells_lt(:,tt)+1)'; % firing probability at any given time for every sarc
    for c = 1 : cell_num
        P(c) = P(c) - prob_increase(flag, tt, distance_matrix, c, cell_num, t0, thres_dict, cell_info);
    end
    
    sarcFire = firingEvent > P;  % sarcomere firing logical vector
    for c = 1 : cell_num
        if sarcFire(c) == 0
            flag(c) = tt;
        end
    end
    cells_lt(:,tt+1) = ((cells_lt(:,tt) + 1) .* sarcFire); % calculate Lt
end
    
    


% probability of firing as a function of time (linear increase)
function [prob, TimeProb] = ProbabilityFunction(t0, t1, tProb, p0)
    TimeProb = 1 : tProb;
    slope = p0/(t1-t0); % slope of how the probability goes up
    prob = zeros(1,tProb); % Probabiliy function
    prob(t1+1 : tProb) = p0; % after time t1, the probability stays at p0
    t_vec = (1 : (t1-t0)); % time vec for slope
    prob(t0+1 : t1) = slope .* t_vec; % from t0 to t1, the probability goes up linearly (relative refractory period)
end

% FN map
function [FN_x_auto, FN_y_auto] = FN_map(FN_x, FN_y, gap, num_segment)
    FN_y_1 = [0,0,FN_y,FN_y]; % y coordinate for every FN segment

    FN_x_auto = zeros(4,num_segment); % initialize x coordinates for all FN segments
    count = 1:num_segment; % from 1 to number of segment
    % creat x coordinates for all FN segments
    FN_x_auto(1,:) = (FN_x+gap)*(count-1);
    FN_x_auto(4,:) = FN_x_auto(1,:);
    FN_x_auto(2,:) = FN_x*count + gap*(count-1);
    FN_x_auto(3,:) = FN_x_auto(2,:);
    % creat y coordinates for all FN segments
    FN_y_auto = repmat(transpose(FN_y_1),1,num_segment);
end

% distance matrix£ºmatrix that contains distance between every two cells
function [distance_matrix] = dist_mat_calculator(cell_center, cell_num)
    distance_matrix = zeros(cell_num,cell_num);
    for i = 1 : cell_num
        for j = 1 : cell_num
            x1 = cell_center(1,i);
            y1 = cell_center(2,i);
            x2 = cell_center(1,j);
            y2 = cell_center(2,j);
            distance_matrix(i,j) = distance_calculator(x1, y1, x2, y2);
        end
    end
end

% cell information
function [cell_info, cell_num, cell_center] = cells(FN_x, FN_y, gap, num_segment, section_w, cell_max_x, cell_min_x)

    n_cell = 1; % initialization for cell number
    num_seg_y_still = 0 : section_w : FN_y; % y coordinated for cells

    for k = 1 : num_segment % loop through every FN segment
        for j = 1 : length(num_seg_y_still)-1 % loop through cell sections
            if rand(1) < 0.5 % left to right
                flag = 0;
                seg_x = 0; % initialize x coordinates of cells in one cell section
                while flag < FN_x % placing a cell on one cell section until the right side of the cell reach the FN map boundary
                    cell_x = rand(1) * (cell_max_x - cell_min_x) + cell_min_x; % length of this cell
                    seg_x = [seg_x,flag + cell_x]; % add the x coordinate of the cell to the seg_x array
                    flag = flag + cell_x; % add cell length to flag
                end    
                seg_x = [seg_x(1:length(seg_x)-1)+(k-1)*(FN_x+gap),(k-1)*(FN_x+gap)+FN_x]; % make the last number of seg_x the FN boundary
            else % from right to left
                flag = FN_x;
                seg_x = FN_x;
                while flag > 0
                    cell_x = rand(1) * (cell_max_x - cell_min_x) + cell_min_x; % length of this cell
                    seg_x = [flag - cell_x,seg_x];
                    flag = flag - cell_x;
                end
                seg_x = [(k-1)*(FN_x+gap),seg_x(2:length(seg_x))+(k-1)*(FN_x+gap)];
            end    
            for i = 1 : length(seg_x)-1  % place cells in one cell section
                cell_info{n_cell} = [seg_x(i),num_seg_y_still(j),seg_x(i+1)-seg_x(i),1000/45,k]; % save the position information of this rectangle into a cell instance
                n_cell = n_cell+1; % go to the next cell
            end
        end 
    end
    cell_num = length(cell_info);
    cell_center = zeros(3,cell_num); % x,y coordinates of cell centers
    for i = 1 : cell_num
        cell_center(1,i) = cell_info{i}(1)+cell_info{i}(3)*0.5;
        cell_center(2,i) = cell_info{i}(2)+cell_info{i}(4)*0.5;
        cell_center(3,i) = cell_info{i}(5);
    end

end

% prob increase caused by cell interaction (linear)
% unit of thres_dist: um
% if distance is longer than thres_dist, there will be no influence
function [prob_in] = prob_increase(flag, current_tt, distance_matrix, cell, cell_num, t0, thres_dict, cell_info)
    prob_in = 0;
    slope = -1/thres_dict;
    if flag(cell) > t0
        for i = 1 : cell_num
            if flag(i) == current_tt-1 && distance_matrix(i,cell) < thres_dict && cell_info{cell}(5) == cell_info{i}(5)
                prob_in = prob_in + distance_matrix(i,cell)*slope + 1;
            end
        end
    end        
end




% distance function
function [distance] = distance_calculator(x1, y1, x2, y2)
    distance = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
end
    