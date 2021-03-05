classdef rt_optimizer < handle
% RT_OPTIMIZER The class for solving optimization tasks
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    properties
        type
        fixed_point_opts = struct(...
            'display', true, ...
            'max_iter', 1e6, ...
            'tol', 1e-8, ...
            'reg_coeff', 0.5 ...
        )
        proximal_ascend_opts = struct(...
            'display', true, ...
            'max_iter', 1e6, ...
            'tol', 1e-8 ...
        )
        auto_rank_opts = struct(...
            'display', true, ...
            'sl', 0.05 ...
        )
    end
    methods
        function obj = rt_optimizer(type)
            switch type
                case {'fixed_point', 'proximal_ascend', 'auto_rank'}
                    obj.type = type;
                otherwise
                    error('RT:OptType', 'Unknown optimizer type: `%s`\n Only `fixed_point`, `proximal_ascend`, and `auto_rank` are available', type);
            end
        end
        function obj = set_options(obj, varargin)
            opts_field = [obj.type, '_opts'];
            for j = 1:2:length(varargin)
                obj.(opts_field).(varargin{j}) = varargin{j+1};
            end
        end
        function [x, info] = run(obj, varargin)
            switch obj.type
                case 'fixed_point'
                    [x, info] = obj.fixed_point(varargin{:});
                case 'proximal_ascend'
                    [x, info] = obj.proximal_ascend(varargin{:});
                case 'auto_rank'
                    [x, info] = obj.auto_rank(varargin{:});
            end
        end
        % ========= Fixed point iteration ============
        function [x, info] = fixed_point(obj, x0, fVal)
            op = obj.fixed_point_opts;
            if op.display
                fprintf('Optimization: fixed point iteration method\n');
                h = rt_fprint('Starting optimization');
            end
            x = x0;
            for j = 1:op.max_iter
                xp = x;
                x = (1-op.reg_coeff) * fVal(x) + op.reg_coeff * xp;
                dx = norm(vec(xp-x));
                stopIter = (dx < op.tol);
                if op.display && (mod(j,op.display) == 0 || j == 1 || stopIter)
                    h = rt_fprint(sprintf('Iteration %d \t\t Delta %.4e', j, dx), h);
                end
                if stopIter
                    break;
                end
            end
            if op.display
                fprintf('\n');
            end
            info.iter = j;
        end
        % ========= Proximal ascend ============
        function [x, info] = proximal_ascend(obj, x0, fFun, fdFun, fProx, LipsConst)
            op = obj.proximal_ascend_opts;
            if op.display
                fprintf('Optimization: proximal gradient ascend\n');
                h = rt_fprint('Starting optimization');
            end
            x1 = x0; z1 = x0; t0 = 0; t1 = 1;
            for j = 1:op.max_iter
                y1 = x1 + t0/t1*(z1-x1) + (t0-1)/t1*(x1-x0);
                z2 = fProx(y1 + fdFun(y1) / LipsConst);
                v2 = fProx(x1 + fdFun(x1) / LipsConst);
                logL_z2 = fFun(z2);
                logL_v2 = fFun(v2);
                if logL_z2 >= logL_v2
                    x2 = z2;
                else
                    x2 = v2;
                end
                x0 = x1; x1 = x2; z1 = z2;
                t0 = t1; t1 = (sqrt(4*t0^2+1)+1)/2;
                dx = norm(vec(x1 - x0));
                stopIter = (dx < op.tol);
                if op.display && (mod(j,op.display) == 0 || j == 1 || stopIter)
                    h = rt_fprint(sprintf('Iteration %d \t\t Delta %.4e', j, dx), h);
                end
                if stopIter
                    break;
                end
            end
            if op.display
                fprintf('\n');
            end
            x = x1;
            info.iter = j;
        end
        % ========= Automatic rank ============
        function [x, info] = auto_rank(obj, rmax, fData)
            op = obj.auto_rank_opts;
            if op.display
                fprintf('=== Automatic rank estimation ===\n');
            end
            pvalRed = false;
            info = deal(cell(1, rmax));
            for r = 1:rmax
                if op.display
                    fprintf('=> Try rank %d\n', r);
                end
                info{r} = fData(r);
                if isnan(info{r}.pval) || info{r}.pval >= op.sl
                    break;
                elseif r > 1 && info{r}.pval < info{r-1}.pval
                    pvalRed = true;
                    r = r - 1;
                    break;
                end
            end
            if op.display
                if info{r}.pval >= op.sl
                    fprintf('=> Rank %d is statistically significant at significance level %s. Procedure terminated.\n', r, num2str(op.sl));
                elseif pvalRed
                    fprintf('=> P-value is maximal (%s) for rank %d. Procedure terminated.\n', num2str(info{r}.pval), r);
                else
                    fprintf('=> Failed to determine optimal rank. Maximal rank %d is taken.\n', rmax);
                end
            end
            x = info{r};
        end
    end
end

