classdef rt_optimizer
    %RT_OPTIMIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        fixed_point_opts = struct(...
            'display', true, ...
            'dispfreq', 50, ...
            'max_iter', 1e6, ...
            'tol', 1e-8, ...
            'reg_coeff', 0.5 ...
        )
        proximal_descend_opts = struct(...
            'display', true, ...
            'dispfreq', 50, ...
            'max_iter', 1e6, ...
            'tol', 1e-8 ...
        )
    end
    
    methods
        function obj = rt_optimizer(type)
            switch type
                case {'fixed_point', 'proximal_descend'}
                    obj.type = type;
                otherwise
                    error('RT:OptType', 'Unknown optimizer type: `%s`\n Only `fixed_point` and `proximal_descend` are available', type);
            end
        end
        
        function obj = set_options(obj, varargin)
            for j = 1:2:length(varargin)
                switch obj.type
                    case 'fixed_point'
                        obj.fixed_point_opts.(varargin{j}) = varargin{j+1};
                    case 'proximal_descend'
                        obj.proximal_descend_opts.(varargin{j}) = varargin{j+1};
                end
            end
        end
        
        function [x, info] = run(obj, x0, varargin)
            switch obj.type
                case 'fixed_point'
                    [x, info] = obj.fixed_point(x0, varargin{:});
                case 'proximal_descend'
                    [x, info] = obj.proximal_descend(x0, varargin{:});
            end
        end
        
        % ========= Fixed point iterations ============
        function [x, info] = fixed_point(obj, x0, fVal)
            op = obj.fixed_point_opts;
            if op.display
                fprintf('Optimization: fixed point iterations method\n');
                h = rt_fprint('Starting optimization');
            end
            x = x0;
            for i = 1:op.max_iter
                xp = x;
                x = (1-op.reg_coeff) * fVal(x) + op.reg_coeff * xp;
                dx = norm(vec(xp-x));
                stopIter = (dx < op.tol);
                if op.display && (mod(i,op.dispfreq) == 0 || i == 1 || stopIter)
                    h = rt_fprint(sprintf('Iteration %d \t\t Delta %.4e', i, dx), h);
                end
                if stopIter
                    break;
                end
            end
            if op.display
                fprintf('\n');
            end
            info.iter = i;
        end
        
        % ========= Proximal descend ============
        function [x, info] = proximal_descend(obj, x0, fLogL, fdLogL, fProximal, LipschitzConstant)
            op = obj.proximal_descend_opts;
            if op.display
                fprintf('Optimization: proximal gradient descend\n');
                h = rt_fprint('Starting optimization');
            end

            x1 = x0; z1 = x0; t0 = 0; t1 = 1; logL_x1 = 0;
            for i = 1:op.max_iter
                y1 = x1 + t0/t1*(z1-x1) + (t0-1)/t1*(x1-x0);
                z2 = fProximal(y1 + fdLogL(y1) / LipschitzConstant);
                v2 = fProximal(x1 + fdLogL(x1) / LipschitzConstant);
                logL_z2 = fLogL(z2);
                logL_v2 = fLogL(v2);
                if logL_z2 >= logL_v2
                    x2 = z2; logL_x2 = logL_z2;
                else
                    x2 = v2; logL_x2 = logL_v2;
                end

                x0 = x1; x1 = x2; z1 = z2;
                t0 = t1; t1 = (sqrt(4*t0^2+1)+1)/2;

                dx = norm(vec(x1 - x0));
                dlogL = logL_x2 - logL_x1;
                logL_x1 = logL_x2;
                stopIter = (dx < op.tol);

                if op.display && (mod(i,op.dispfreq) == 0 || i == 1 || stopIter)
                    h = rt_fprint(sprintf('Iteration %d \t\t Delta %.4e \t\t dLogL %.4e', i, dx, dlogL), h);
                end
                if stopIter
                    break;
                end
            end
            if op.display
                fprintf('\n');
            end
            x = x1;
            info.iter = i;
        end
    end
end

