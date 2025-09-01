% Example of using binary_reader_wrapper
% Need to change "format" to 2 when post-processing the example cases

clear; clc; close all;

mfcPath = '/disk/simulations/PhaseChange/ShockTube/1D/StrongCollapse/6Eqn/';

%% Fluid properties

% Water Liquid (1) and air(2)
pii = [1E9 0]'; % p_infty, pa
q = [-1.167E6 0E6]'; % J/Kg
qp = [0 0]'; % J/KgK
cv = [1816 717.5]'; % J/KgK
cp = [4267 1006]'; % J/KgK
gama = cp ./ cv ;

%% 1D data extraction
RelMod = {'p', 'pT'} ;
DiscLevel = {'N1E3', 'N2E3'} ;

for rm = 1:length(RelMod)
    for dl = 1:length(DiscLevel)

        binDir = fullfile(mfcPath, RelMod{rm}, DiscLevel{dl}, 'binary' ) ;
        
        [alpha_rho1, alpha_rho2, mom1, vel1, E, alpha_rho_e1, alpha_rho_e2, pres, tCoord, xCoord] = binary_reader_wrapper(binDir, 1) ;
        
        tCoord = repmat( tCoord, length( xCoord ) , 1 ) ;
        
        %% calculations for the x-t plots
        
        alpha1 = (gama(1) - 1) .* (alpha_rho_e1 - alpha_rho1 .* q(1) ) ./ ( pres + gama(1) .* pii(1) ) ;
        alpha2 = (gama(2) - 1) .* (alpha_rho_e2 - alpha_rho2 .* q(2) ) ./ ( pres + gama(2) .* pii(2) ) ;
        
        rho1 = alpha_rho1 ./ alpha1 ;
        rho2 = alpha_rho2 ./ alpha2 ;
        
        T1 = ( pres + pii(1) ) ./ ( ( gama(1) - 1 ) .* cv(1) .* rho1 ) ;
        T2 = ( pres + pii(2) ) ./ ( ( gama(2) - 1 ) .* cv(2) .* rho2 ) ;
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        tl = tiledlayout( 2, 2, 'TileSpacing', 'compact' );
        fs = 30 ;
        
        Lx = xCoord(end,1) - xCoord(1,1) ;        

        nexttile(1) ;
        contourf(tCoord ./ tCoord(end), xCoord ./ Lx, pres) ;
        title( '$ p \; [Pa] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.1f');
        ax = gca; ax.FontSize = fs;
        
        nexttile(2) ;
        contourf(tCoord ./ tCoord(end), xCoord ./ Lx, T2) ;
        title( '$ T_{v} \; [K] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.1f');
        ax = gca; ax.FontSize = fs;

        nexttile(3) ;
        contourf(tCoord ./ tCoord(end), xCoord ./ Lx, E) ;
        title( '$ E \; [J] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.1f');
        ax = gca; ax.FontSize = fs;

        nexttile(4) ;
        contourf(tCoord ./ tCoord(end), xCoord ./ Lx, alpha_rho2) ;
        title( '$ m_{2} \; [kg/m^3] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.1f');
        ax = gca; ax.FontSize = fs;

        title( tl, strcat(RelMod{rm}, '-relax. Discr ', DiscLevel{dl}), 'interpreter', 'latex', 'Fontsize', fs);
        xlabel( tl, '$ t \; [s] $', 'interpreter', 'latex', 'Fontsize', fs);
        ylabel( tl, '$x/L_x$', 'Interpreter', 'latex', 'FontSize', fs ) ;
        
        % colormap( flipud(gray(256)) ) ;
    
        clearvars alpha_rho1 alpha_rho2 mom1 vel1 E alpha_rho_e1 alpha_rho_e2 pres tCoord xCoord
        
        savefig(fig, fullfile('/disk/simulations/PhaseChange/ShockTube/1D/StrongCollapse/6Eqn/Figures/', strcat(RelMod{rm}, DiscLevel{dl} ) ), '-v7.3' );

        close 
    end
end