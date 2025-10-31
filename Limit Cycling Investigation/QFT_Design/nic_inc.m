function nic_inc
% nic_inc
% To draw an inverse nichols chart, you need to flip the nichols chart in 
% the installed nicchart.m. The code below does this for you by clearing
% the function and using N or I for selecting which chart is required. You
% have to provde a directory in the path in which the code files are
% stored. After that, functions like grid, ngrid and nichols will use the
% chart you have selected

% Change this directory to suite your set-up (i.e. edit this line!)
NicDirectory ='C:\Users\Zwivh\OneDrive - University of Cape Town\5th Year\2nd Semester\EEE4022S\MATLAB\Controller Design (Water Level)\QFT_Design';  

n = input('N for Nichols, I for INC ','s');
n=upper(n)
WorkingDirectory=cd;          % find current directory
cd (NicDirectory);          % go to where the files are
clear nicchart      % clear the internal version of the function
switch n
    case 'N'
        disp('Nichols chart selected')
        [status,message] = copyfile('nicchart_src.m','nicchart.m');
    case 'I'
        disp('Inverse Nichols chart selected')
        [status,message] = copyfile('innicchart_src.m','nicchart.m');
     otherwise
        disp('other value')
end
cd(WorkingDirectory);   % go back to where you were

end

