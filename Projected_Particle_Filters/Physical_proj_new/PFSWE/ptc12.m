
function distinctColors = ptc12(n, varargin)
% PTC12 color palettes from Paul Tol at SRON
% Creates a palette of n (between 1 and 12) distinctive
% colors that are:
%  - distinct for all people, including colour-blind readers
%  - distinct from black and white
%  - distinct on screen and paper
%  - still match well together
%
% Output is (n,3) RGB colors
% A light grey is included at n+1 (can indicate absence of data)
%
% Using n=0 gives an alternative 4 color palette 
% that is optimized for printing in grey scale 
% 
% calling ptc(n,'check') will show a figure displaying 
% the color scheme that is selected
% Colors from ptc6 are slightly brighter than those from ptc12
%
% All colors and schemes created by Paul Tol at SRON:
% https://personal.sron.nl/~pault/
% Written by Nans Bujan, PhD, 02/03/2019. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n <0 || n >12
    error('N must be comprised between 0 and 12')
end    

% Colors definition in RGB
PAULTOLCOLORS =     [51,    34,     136;...                  % 1
                    102,    153,    204;...
                    136,    204,    238;...
                    68,     170,    153;...
                    17,     119,    51;...                   % 5
                    153,    153,    51;...
                    221,    204,    119;...
                    102,    17,     0;...
                    204,    102,    119;...
                    170,    68,     102;...                  % 10
                    136,    34,     85;...
                    170,    68,     153;...
                    68,     119,    170;...
                    221,    221,    221]...                  % light grey
                    / 255; 

% Colors schemes 
PAULTOLSCHEMES = {[2, 0, 5, 11], ... % 4 color scheme for grey scale % 1
                [12], ...                                       
                [12, 8], ... 
                [12, 6, 8], ... 
                [12, 4, 6, 8], ...                           % 5 
                [0, 2, 4, 6, 8], ...                             
                [0, 2, 4, 6, 8, 11], ... 
                [0, 2, 3, 4, 6, 8, 11], ... 
                [0, 2, 3, 4, 5, 6, 8, 11], ... 
                [0, 2, 3, 4, 5, 6, 8, 10, 11], ...           % 10
                [0, 2, 3, 4, 5, 6, 1, 8, 10, 11], ...             
                [0, 7, 2, 3, 4, 5, 6, 1, 8, 10, 11], ...
                [0, 7, 2, 3, 4, 5, 6, 1, 8, 9, 10, 11]};

distinctColors = PAULTOLCOLORS([PAULTOLSCHEMES{n+1}+1,14],:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays a figure to show the selected colors:
    if strcmp(varargin,'check')
        n(n==0) = 4; % still show 4 colors with the special 4 color scheme
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