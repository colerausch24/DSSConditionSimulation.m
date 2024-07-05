n = 10;
m = round(10.5 * n);
expirements = 1000000000; %1,000,000,000
maxRatio_CprimeW = 0;
minRatio_CprimeW = inf;
sumRatio_CprimeW = 0;
count=0;
ratios_CprimeW = zeros(1, expirements);
ReLU = @(x) max(0,x);
noDSS=0;

for i = 1:expirements
    while true
        try
            W = randn(m, n);
            x = randn(n, 1);
            S = find(W * x >= 0);
            if isempty(S)
                error();
            end
            W_S = W(S, :);
            break;
        catch ME
        end
    end

    sigma_min = min(svd(W_S));
    C_prime_W = sigma_min / sqrt(2);

    x0 = randn(n, 1);
    x1 = randn(n, 1);
    LHS = norm(ReLU(W * x0) - ReLU(W * x1), 2);
    RHS_prime = C_prime_W * norm(x0 - x1, 2);

    ratio_CprimeW = LHS / RHS_prime;
    ratios_CprimeW(i) = ratio_CprimeW;
    maxRatio_CprimeW = max(maxRatio_CprimeW, ratio_CprimeW);
    minRatio_CprimeW = min(minRatio_CprimeW, ratio_CprimeW);
    sumRatio_CprimeW = sumRatio_CprimeW + ratio_CprimeW;

    if ratio_CprimeW <= 1
        count = count + 1;
        fprintf('Need to test DSS condition for this W..Iteration: %d:\n', i);
        [~, noDSS] = test_find_intersections(W, x0, x1, noDSS);
    end
end

fprintf('expirements:  %d\n',expirements);
fprintf('Maximum LHS/RHS_prime Ratio: %f\n', maxRatio_CprimeW);
fprintf('Minimum LHS/RHS_prime Ratio: %f\n', minRatio_CprimeW);
fprintf('Average LHS/RHS_prime Ratio: %f\n', avgRatio_CprimeW);
fprintf('Times LHS/RHS < 1:  %d\n',count);
percent_count=(count/expirements)*100;
fprintf("Percent of times LHS/RHS <= 1: %f\n",percent_count);
percent_count2=1-percent_count;
fprintf("Percent of times LHS/RHS > 1: %f\n",percent_count2);
if noDSS ~= 0
    noDSSpercent = (count / noDSS) * 100;
    fprintf("Percentage of instances where W did not satisfy the inequality and also did not meet the DSS condition: %d\n", noDSSpercent);
else
    fprintf('No cases where W did not satisfy the inequality.\n');
end


% Find boundary points for wedges
function ts = find_intersections(W, x1, x0)
%W is m x n
%x0 and x1 are n x 1
%therefore W*x= (mxn)x(nx1)=(mx1)

    numerator = (-1*W)*x0;
    denominator = W*(x1-x0);
    ts = numerator ./ denominator;
    % Handling error cases
    ts(denominator == 0 | ts < 0 | ts > 1) = NaN;
end

function [passDSS, noDSS] = test_find_intersections(W, x0, x1, noDSS)
    intersection_ts = find_intersections(W, x0, x1);
    
    % Filter out NaN values and ensure ts are within [0, 1]
    intersection_ts = intersection_ts(~isnan(intersection_ts) & intersection_ts >= 0 & intersection_ts <= 1);
    
    % Sort the ts array and include 0 and 1 as boundaries
    intersection_ts = sort([0, reshape(intersection_ts, 1, []), 1]);
    passDSS = true;
    
    % Loop through each segment defined by the ts
    for j = 1:length(intersection_ts) - 1
        % Calculate midpoint
        t = (intersection_ts(j) + intersection_ts(j + 1)) / 2;
        
        % Line segment between x0 and x1 at midpoint t
        l_t = (1 - t) * x0 + (t * x1);
        
        % Check the DSS condition for the midpoint
        if ~has_DSS(W, l_t)
            fprintf('Segment [%f, %f] does not satisfy DSS condition.\n', intersection_ts(j), intersection_ts(j+1));
            passDSS = false;
            noDSS = noDSS + 1;
            break;
        else
            fprintf('Segment [%f, %f] satisfies DSS condition.\n', intersection_ts(j), intersection_ts(j+1));
        end
    end
end


function result = has_DSS(W, l_t)
    S = W * l_t >= 0;
    result = all(S);
end
