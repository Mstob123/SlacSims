% LIAR setting for the FODO beamline

mat_liar('  reset, all');
mat_liar('  set_control, debug   = 0');    % Debug-level               
%mat_liar('               outlun  = 6,');    % Output unit (6 = terminal)
%mat_liar('               outfile = ' '');   % Optional output file name
%
mat_liar('  read_xsif, ');
mat_liar('  file = ''fodo100.mad'',');
mat_liar('             energy = 500.,');
mat_liar('             echo   = .f.,');
mat_liar('             line   = ''fodo'' ');
%
mat_liar('  calc_twiss, ');
mat_liar('              betax   =  20.2264,');   % Initial beta_x  
mat_liar('              betay   =  20.2264,');   % Initial beta_y  
mat_liar('              alphax  = -1.16246,');   % Initial alpha_x 
mat_liar('              alphay  =  1.16246 ');   % Initial alpha_y 
%
mat_liar('  set_initial, x       = 0.d-6,');   % Initial x offset in m      
mat_liar('               angx    = 0.d-6,');   % Initial x beam angle in rad
mat_liar('               y       = 0.-6,');   % Initial y offset in m      
mat_liar('               angy    = 0.d-6,');   % Initial y beam angle in rad
mat_liar('               energy  = 500., ');    % Initial beam energy in GeV 
mat_liar('               espread = 0.5,  '); % Initial energy spread in GeV
mat_liar('               nemitx  = 360d-6,'); % Init. x emit. in rad*m
mat_liar('               nemity  = 4.0d-8,'); % Init. y emit. in rad*m
mat_liar('               betax   =  20.2264,');   % Initial beta_x  
mat_liar('               betay   =  20.2264,');   % Initial beta_y  
mat_liar('               alphax  = -1.16246,');   % Initial alpha_x 
mat_liar('               alphay  =  1.16246 ');   % Initial alpha_y 
%
mat_liar('  set_beam, current = 7.5d9,'); % Number of part. per bunch
mat_liar('            blength = 100.d-6,');  % RMS bunch length in m 
mat_liar('            ecut    =  3,     ');  % Energy cut in sigmas
mat_liar('            zcut    =  3,     ');  % Longit. cut in sigmas
mat_liar('            nb      = 1,      ');  % Nb. of bunches
mat_liar('            ns      = 10,     ');  % Nb. of slices
mat_liar('            nm      = 5,      ');  % Nb. of macro-particles
mat_liar('            bspace  = 4.2     ');  % Spacing in m between bunches
% 
% Definition of the girder supports:
%mat_liar('define_support,ngirder=0, qsupport= .t. ');
mat_liar('define_support,ngirder=1 ');
%
%mat_liar('show_definitions  ');
%mat_liar('show_misalign  ');
%mat_liar('show_support  ');
%
return
