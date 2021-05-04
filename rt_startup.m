% RT_STARTUP Includes the library directories in the search paths. Run before using the library.
% Documentation: https://github.com/PQCLab/mRootTomography/blob/master/Documentation.md
% The code is licensed under GPL v3
% Author: Boris Bantysh, 2021
addpath(pwd, strcat(pwd, '\Tools'));
addpath(pwd, strcat(pwd, '\Statistics'));
fprintf('All paths are set. Library is ready to use.\n');
fprintf('For documentation see library page on GitHub: <a href="https://github.com/PQCLab/mRootTomography">PQCLab/mRootTomography</a>.\n');
fprintf('Some examples could be found in Examples directory.\n');