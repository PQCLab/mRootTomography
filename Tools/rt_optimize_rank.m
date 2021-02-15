function [data, r, data_r] = rt_optimize_rank(dim, sl, display, fData)
%RT_OPTIMIZE_RANK Summary of this function goes here
%   Detailed explanation goes here

if display
    fprintf('=== Automatic rank estimation ===\n');
end

pvalRed = false;
data_r = deal(cell(1,dim));
for r = 1:dim
    if display
        fprintf('=> Try rank %d\n', r);
    end
    data_r{r} = fData(r);
    if isnan(data_r{r}.pval) || data_r{r}.pval >= sl
        break;
    elseif r > 1 && data_r{r}.pval < data_r{r-1}.pval
        pvalRed = true;
        r = r - 1;
        break;
    end
end

if display
    if data_r{r}.pval >= sl
        fprintf('=> Rank %d is statistically significant at significance level %s. Procedure terminated.\n', r, num2str(sl));
    elseif pvalRed
        fprintf('=> P-value is maximal (%s) for rank %d. Procedure terminated.\n', num2str(data_r{r}.pval), r);
    else
        fprintf('=> Failed to determine optimal rank. Maximal rank %d is taken.\n', r);
    end
end

data = data_r{r};

end

