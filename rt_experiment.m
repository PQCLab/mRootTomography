classdef rt_experiment < matlab.mixin.Copyable
% RT_EXPERIMENT The class for working with the quantum tomography data
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
    properties
        dim                 % Hilbert space dimension
        obj_type            % Object type (`'state'` or `'process'`)
        stat_pkg = 'auto'   % Statistics package
        proto = {}          % Measurements operators
        nshots = []         % Measurements repetitions
        clicks = {}         % Number of observed measurements outcomes
        vec_proto = []      % Matrix form of the whole measurements protocol
        vec_nshots = []     % Matrix form of the measurements repetitions
        vec_clicks = []     % Matrix form of the number of observed measurements outcomes
    end
    methods
        function obj = rt_experiment(dim, obj_type, stat)
            obj.dim = dim;
            if nargin > 2
                if ischar(stat)
                    stat = lower(stat);
                    st_allowed = [reshape(fieldnames(rt_stat.buildin), 1, []), 'auto'];
                    if ~any(strcmp(stat, st_allowed))
                        error('RT:StatsType', 'Unknown statistics type: `%s`\n Available types: %s', stat, strjoin(st_allowed, ', '));
                    end
                end
                obj.stat_pkg = stat;
            end
            switch obj_type
                case {'state', 'process'}
                    obj.obj_type = lower(obj_type);
                otherwise
                    error('RT:ObjectType', 'Unknown object type: `%s`\n Only `state`, `process` are available', obj_type);
            end
        end
        function obj = set_data(obj, varargin)
            for j = 1:2:length(varargin)
                field = lower(varargin{j});
                obj.(field) = varargin{j+1};
                if strcmp(field, 'proto')
                    % Check whether the process protocol is set as
                    % prepare-measure pairs and transform to chi-matrix
                    % measurement operators
                    if strcmp(obj.obj_type, 'process') && iscell(obj.proto) && all(size(obj.proto) > 1)
                        if size(obj.proto, 2) == 2
                            obj.proto = transpose(obj.proto);
                        end
                        obj.proto = cellfun(@(prep, meas) rt_kron3d(conj(prep), meas), obj.proto(1,:), obj.proto(2,:), 'UniformOutput', false);
                    end
                    % Check whether the protocol is set as 3d array and
                    % transform it ti cell array
                    if ~iscell(obj.proto)
                        dims = size(obj.proto);
                        obj.proto = mat2cell(obj.proto, dims(1), dims(2), ones(1, dims(3)));
                    end
                    obj.proto = reshape(obj.proto, 1, []);
                    obj.vec_proto = [];
                    continue;
                end
                if strcmp(field, 'clicks')
                    if ~iscell(obj.clicks)
                        obj.clicks = reshape(num2cell(obj.clicks), 1, []);
                    end
                    obj.vec_clicks = [];
                    continue;
                end
                if strcmp(field, 'nshots')
                    obj.vec_nshots = [];
                end
            end
            if ~isempty(obj.nshots) && ~isempty(obj.proto)
                if length(obj.nshots) == 1
                    obj.nshots = rt_nshots_divide(obj.nshots, size(obj.proto, 2));
                elseif length(obj.nshots) ~= size(obj.proto, 2)
                    error('RT:ExpNumberMismatch', 'Length of nshots array does not match length of proto array');
                end
            end
        end
        function data = get_field(obj, field)
            if strcmp(field, 'vec_proto') && isempty(obj.vec_proto)
                obj.vec_proto = rt_meas_matrix(cat(3, obj.proto{:}));
            elseif strcmp(field, 'vec_nshots') && isempty(obj.vec_nshots)
                n = arrayfun(@(nj,Mj) ones(size(Mj{1},3),1)*nj, obj.nshots(:), obj.proto(:), 'UniformOutput', false);
                obj.vec_nshots = vertcat(n{:});
            elseif strcmp(field, 'vec_clicks') && isempty(obj.vec_clicks)
                k = cellfun(@(kj) kj(:), obj.clicks, 'UniformOutput', false);
                obj.vec_clicks = vertcat(k{:});
            end
            data = obj.(field);
        end
        function stat = stat(obj)
            if ischar(obj.stat_pkg)
                if strcmp(obj.stat_pkg, 'auto')
                    if any(isinf(obj.nshots))
                        obj.nshots = ones(size(obj.nshots));
                        obj.stat_pkg = 'asymptotic';
                    else
                        imat = eye(obj.dim);
                        if all(cellfun(@(pr) size(pr, 3) == 1, obj.proto)) % is binomial
                            obj.stat_pkg = 'binomial';
                        elseif strcmp(obj.obj_type, 'state') && all(cellfun(@(pr) norm(sum(pr, 3) - imat) < 1e-5, obj.proto)) % is state povm
                            obj.stat_pkg = 'polynomial';
                        elseif strcmp(obj.obj_type, 'process') && all(cellfun(@(pr) norm(rt_prttrace(sum(pr, 3), [obj.dim, obj.dim], 1) - imat) < 1e-5, obj.proto)) % is process povm
                            obj.stat_pkg = 'polynomial';
                        else
                            error('RT:StatsTypeAuto', 'Failed to determine statistics type. Please, specify `stat` manually.');
                        end
                    end
                end
                obj.stat_pkg = rt_statistics.get(obj.stat_pkg);
            end
            stat = obj.stat_pkg;
        end
        function p = get_probs_dm(obj, dm, tol)
            if nargin < 3
                tol = 0;
            end
            p = abs(obj.get_field('vec_proto') * dm(:));
            p(p < tol) = tol;
        end
        function p = get_probs_sq(obj, sq, tol)
            if nargin < 3
                tol = 0;
            end
            p = obj.get_probs_dm(sq * sq', tol);
        end
        function [n, k, p] = nkp(obj, dm)
            n = obj.get_field('vec_nshots');
            k = obj.get_field('vec_clicks');
            if nargin > 1
                p = obj.get_probs_dm(dm, 1e-15);
            end
        end
        % ========= Sampling ============
        function obj = simulate(obj, dm)
            clk = cell(size(obj.proto));
            for j = 1:numel(obj.proto)
                probs = abs(rt_meas_matrix(obj.proto{j}) * dm(:));
                clk{j} = obj.stat().sample(obj.nshots(j), probs);
            end
            obj.vec_clicks = [];
            obj.set_data('clicks', clk);
        end
        % ========= Likelihood ============
        function f = logL_dm(obj, dm)
            [n, k, p] = obj.nkp(dm);
            f = obj.stat().logL(n, k, p);
        end
        function f = logL_sq(obj, sq)
            f = obj.logL_dm(sq * sq');
        end
        function df = dlogL_sq(obj, sq)
            [n, k, p] = obj.nkp(sq * sq');
            b = obj.stat().dlogL(n, k, p);
            B = obj.get_field('vec_proto');
            df = 2 * reshape(B' * b, size(sq, 1), []) * sq;
        end
        function mu = logL_eq_mu(obj)
            [n, k] = obj.nkp();
            mu = obj.stat().logL_mu(n, k);
        end
        function jmat = logL_eq_jmat_dm(obj, dm)
            [n, k, p] = obj.nkp(dm);
            [b, b0] = obj.stat().logL_jmat(n, k, p);
            B = obj.get_field('vec_proto');
            jmat = reshape(B' * b, obj.dim, obj.dim);
            if b0 ~= 0
                jmat = jmat + b0 * eye(obj.dim);
            end
        end
        % ========= Chi-squared ============
        function chi2 = chi2_dm(obj, dm)
            [n, k, p] = obj.nkp(dm);
            chi2 = obj.stat().chi2(n, k, p);
        end
        function nu = deg_f_rank(obj, obj_rank)
            nu = obj.stat().deg_f(obj.get_field('clicks'));
            if strcmpi(obj.obj_type, 'state')
                nu_dm = (2*obj.dim - obj_rank)*obj_rank - 1;
            elseif strcmpi(obj.obj_type, 'process')
                dim2 = obj.dim^2;
                nu_dm = (2*dim2 - obj_rank)*obj_rank - dim2;
            end
            nu = nu - nu_dm;
        end
    end
end

