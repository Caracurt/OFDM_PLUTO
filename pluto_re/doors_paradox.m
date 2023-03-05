% paradox of 3 doors

num_el = 3;
exp_num = 10000;

% EXP1 // success
prob_success_no_change = [];
prob_success_change = [];

for exp_idx = 1:exp_num

    doors = zeros(num_el, 1);
    pos_gold = randi([1, num_el], 1, 1);
    doors(pos_gold) = 1;

    % init choice random
    c1 = randi([1, num_el], 1, 1);

    % position exept choice
    pos_ex = [1:c1-1, c1+1:num_el];

    % check positions
    pos_reveal = find(doors(pos_ex) == 0);
    pos_reveal = pos_ex(pos_reveal(1));

    assert(doors(pos_reveal) == 0, 'Error in Reveal');
    assert(pos_reveal ~= c1, 'Error in reveal');

    c2_in = find(pos_ex ~= pos_reveal);
    c2 = pos_ex(c2_in);

    assert(c1 ~= c2, 'Error in final choice');

    
    prob_success_no_change = [prob_success_no_change, double(c1 == pos_gold)];
    prob_success_change = [prob_success_change, double(c2 == pos_gold)];
end

prob_success_no_change = sum(prob_success_no_change) / length(prob_success_no_change);
prob_success_change = sum(prob_success_change) / length(prob_success_change);

fprintf('EXP1 ProbWin: NoChangeStrat=%f ProbWin: ChangeStrat=%f\n', prob_success_no_change, prob_success_change);


% EXP2 // new knowledge
prob_success_no_change = [];
prob_success_change = [];

for exp_idx = 1:exp_num

    doors = zeros(num_el, 1);
    pos_gold = randi([1, num_el], 1, 1);
    doors(pos_gold) = 1;

    % init choice random
    c1 = randi([1, num_el], 1, 1);

    % position exept choice
    pos_ex = [1:c1-1, c1+1:num_el];

    % check positions
    % random reveal
    pos_reveal_int = randi([1, num_el-1], 1, 1);

    pos_reveal = pos_ex(pos_reveal_int);

    if doors(pos_reveal) == 0
        c2_in = find(pos_ex ~= pos_reveal);
        c2 = pos_ex(c2_in);

        assert(c1 ~= c2, 'Error in final choice');
    else
        c2 = c1;
    end

    %pos_reveal = find(doors(pos_ex) == 0);
    %pos_reveal = pos_ex(pos_reveal(1));

    %assert(doors(pos_reveal) == 0, 'Error in Reveal');
    %assert(pos_reveal ~= c1, 'Error in reveal');

    %c2_in = find(pos_ex ~= pos_reveal);
    %c2 = pos_ex(c2_in);

    %

    
    prob_success_no_change = [prob_success_no_change, double(c1 == pos_gold)];
    prob_success_change = [prob_success_change, double(c2 == pos_gold)];
end

prob_success_no_change = sum(prob_success_no_change) / length(prob_success_no_change);
prob_success_change = sum(prob_success_change) / length(prob_success_change);

fprintf('EXP2 ProbWin: NoChangeStrat=%f ProbWin: ChangeStrat=%f\n', prob_success_no_change, prob_success_change);

