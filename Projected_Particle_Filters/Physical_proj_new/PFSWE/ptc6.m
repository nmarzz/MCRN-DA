
function distinctColors = ptc6(n, varargin)
% PTC6: 'Bright' color palettes from Paul Tol at SRON
% Creates a palette of n (between 1 and 6) distinctive
% colors that are:
%  - distinct for all people, including colour-blind readers
%  - distinct from black and white
%  - distinct on screen and paper
%  - still match well together
%
% Output is (n,3) RGB colors
% A light grey is included at n+1 (can indicate absence of data)
% 
% Calling ptc6(n,'check') will show a figure displaying 
% the color scheme that is selected
% Colors from ptc6 are slightly brighter than those from ptc12
%
% All colors and schemes created by Paul Tol at SRON:
% https://personal.sron.nl/~pault/
% Written by Nans Bujan, PhD, 02/03/2019. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n <0 || n >6
    error('N must be comprised between 0 and 6')
end    

% Colors definition in RGB
PAULTOLCOLORS =     [68,    119,    170;...                 % blue
                    102,    204,    238;...                 % cyan   
                     34,    136,     51;...                 % green
                    204,    187,     68;...                 % yellow
                    238,    102,    119;...                 % red
                    170,     51,    119;...                 % pink
                    187,    187,    187]...                 % light grey
                    / 255; 

% Colors schemes 
PAULTOLSCHEMES = { ...
                [1], ...                                       
                [1, 5], ... 
                [1, 5, 3], ... 
                [1, 5, 3, 4], ...                           
                [1, 2, 3, 4, 5], ...                             
                [1, 2, 3, 4, 5, 6]};

distinctColors = PAULTOLCOLORS([PAULTOLSCHEMES{n},7],:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays a figure to show the selected colors:
    if strcmp(varargin,'check')
        X       = [0:0.1:5]; % a vector to plot some SIN functions
        LW      = 6; % Width of plot lines
        figure
        
        sp1     = subplot(2,1,1);   % white background
            for  i = 1:n+1
                plot(X,sin(X+i/2),'Color',distinctColors(i,:),'LineWidth',LW)
                hold on
            end
            grid on
            title(['Showing ',num2str(n),' colors (and one light grey)'])
        
        sp2     = subplot(2,1,2);   % black background
            copyobj(get(sp1,'Children'),sp2);
            set(gca, 'color', 'k')
            grid on
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end