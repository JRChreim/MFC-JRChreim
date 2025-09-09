% Example of using binary_reader_wrapper
% Need to change "format" to 2 when post-processing the example cases

clear; clc; close all;

% toggle location
loc = 'carpenter';
switch loc
    case 'carpenter'
        mfcPath = '/p/global/jrchreim/simulations/PhaseChange/1D/BubbleDynamics/StrongCollapse/6Eqn/';
        RelMod = {'pTFinalTest', 'pFinalTest'} ;
        % RelMod = {'pTFinalTest'} ;
        DiscLevel = {'N160E3', 'N320E3', 'N640E3', 'N1280E3'} ;
        % DiscLevel = {'N160E3'} ;
        compliment = 'Cartesian/BC-6/C000E-00';
        FigFolder = '/p/global/jrchreim/Figures/';
    case 'local'
        mfcPath = '/disk/simulations/PhaseChange/ShockTube/1D/StrongCollapse/6Eqn/';
        RelMod = {'pT', 'p'} ;
        DiscLevel = {'N1E3', 'N2E3'} ;
        compliment = '';
        FigFolder = '/disk/simulations/PhaseChange/ShockTube/1D/StrongCollapse/6Eqn/Figures/';
end

%% Fluid properties

% Water Liquid (1) and air(2)
pii = [1E9 0]'; % p_infty, pa
q = [-1.167E6 0E6]'; % J/Kg
qp = [0 0]'; % J/KgK
cv = [1816 717.5]'; % J/KgK
cp = [4267 1006]'; % J/KgK
gama = cp ./ cv ;


for rm = 1:length(RelMod)
    for dl = 1:length(DiscLevel)

        binDir = fullfile(mfcPath, RelMod{rm}, compliment, DiscLevel{dl}, 'binary' ) ;
        
        [alpha_rho1, alpha_rho2, mom1, vel1, E, alpha_rho_e1, alpha_rho_e2, pres, tCoord, xCoord] = binary_reader_wrapper(binDir, 1) ;
        
        tCoord = repmat( tCoord, length( xCoord ) , 1 ) ;
        
        %% calculations for the x-t plots
        
        alpha1 = (gama(1) - 1) .* (alpha_rho_e1 - alpha_rho1 .* q(1)) ./ ( pres + gama(1) .* pii(1) ) ;
        alpha2 = (gama(2) - 1) .* (alpha_rho_e2 - alpha_rho2 .* q(2)) ./ ( pres + gama(2) .* pii(2) ) ;
        
        rho1 = alpha_rho1 ./ alpha1 ;
        rho2 = alpha_rho2 ./ alpha2 ;
        
        T1 = ( pres + pii(1) ) ./ ( ( gama(1) - 1 ) .* cv(1) .* rho1 ) ;
        T2 = ( pres + pii(2) ) ./ ( ( gama(2) - 1 ) .* cv(2) .* rho2 ) ;
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        tl = tiledlayout( 2, 2, 'TileSpacing', 'compact' );
        fs = 30 ;
   
        % SETTING UP minimum and maximum values for the colorbars
        if ( dl == 1 && rm == 1)
            maxP = max( max( pres ) ) ; minP = max( min( pres ) ) ;
            maxT2 = max( max( T2 ) )  ; minT2 = max( min( T2 ) ) ;
            maxE = max( max( E ) ) ; minE = max( min( E ) ) ;
            maxM2 = max( max( alpha_rho2 ) ) ; minM2 = max( min( alpha_rho2 ) ) ;
        end

        Lx = xCoord(end,1) - xCoord(1,1) ;        

        tOtend = tCoord ./ tCoord(end) ;
        xOL = xCoord ./ Lx ;

        xC{1} = (tOtend > 0 & tOtend < 0.25) ;
        % xC{1} = (xOL >= 0.0 & xOL <= 1.0) ;

        nexttile(1) ;
        contourf(reshape(tOtend(xC{1}), size(xOL, 1), []), reshape(xOL( xC{1} ), size(xOL, 1), []), reshape(pres(xC{1}), size(xOL,1), [] ) ) ;
        % contourf(reshape(tOtend(xC{1}), [], size(xOL, 2)), reshape(xOL( xC{1} ), [], size(xOL, 2)), reshape(pres(xC{1}), [], size(xOL,2) ) ) ;
        title( '$ p \; [Pa] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.2f'); ytickformat('%.2f');
        text(0.1,0.6,['p_{max}: ' num2str( max( max( pres ) ) ) ], 'FontSize', fs);
        text(0.1,0.4,['p_{min}: ' num2str( min( min( pres ) ) ) ], 'FontSize', fs);
        ax = gca; ax.FontSize = fs;

        clim( [minP, maxP] ) ;
        colorbar ; set(gca,'ColorScale','log') ;

        % xC{2} = (xOL > 0.2 & xOL < 0.7) ;
        xC{2} = (xOL >= 0.0 & xOL <= 1.0) ;
        nexttile(2) ;
        contourf(reshape(tOtend(xC{2}), [], size(xOL, 2)), reshape(xOL( xC{2} ), [], size(xOL, 2)), reshape(T2(xC{2}), [], size(xOL,2) ) ) ;
        title( '$ T_{v} \; [K] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.2f');
        text(0.5,0.6,['T_{max}: ' num2str( max( max( T2 ) ) ) ], 'FontSize', fs);
        text(0.5,0.4,['T_{min}: ' num2str( min( min( T2 ) ) ) ], 'FontSize', fs);
        ax = gca; ax.FontSize = fs;

        clim( [minT2, maxT2] ) ;
        colorbar ;

        xC{3} = (xOL > 0.32 & xOL < 0.34) ;
        % xC{3} = (xOL >= 0.0 & xOL <= 1.0) ;
        nexttile(3) ;
        contourf(reshape(tOtend(xC{3}), [], size(xOL, 2)), reshape(xOL( xC{3} ), [], size(xOL, 2)), reshape(E(xC{3}), [], size(xOL,2) ) ) ;
        title( '$ E \; [J] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.4f');
        text(0.5,0.335,['E_{max}: ' num2str( max( max( E ) ) ) ], 'FontSize', fs);
        text(0.5,0.325,['E_{min}: ' num2str( min( min( E ) ) ) ], 'FontSize', fs);
        ax = gca; ax.FontSize = fs;

        clim( [minE, maxE] ) ;
        colorbar ; set(gca,'ColorScale','log') ;

        xC{4} = (xOL > 0.32 & xOL < 0.34) ;
        % xC{4} = (xOL >= 0.0 & xOL <= 1.0) ;
        nexttile(4) ;        
        contourf(reshape(tOtend(xC{4}), [], size(xOL, 2)), reshape(xOL( xC{4} ), [], size(xOL, 2)), reshape(alpha_rho2(xC{4}), [], size(xOL,2) ) ) ;
        title( '$ m_{2} \; [kg/m^3] $', 'interpreter', 'latex', 'Fontsize', fs);
        text(0.5,0.325,['m_{2,max}: ' num2str( max( max( alpha_rho2 ) ) ) ], 'FontSize', fs);
        text(0.5,0.335,['m_{2,min}: ' num2str( min( min( alpha_rho2 ) ) ) ], 'FontSize', fs);
        xtickformat('%.1f'); ytickformat('%.4f');
        ax = gca; ax.FontSize = fs;

        clim( [minM2, maxM2] ) ;
        colorbar ; set(gca,'ColorScale','log') ;

        title( tl, strcat(RelMod{rm}, '-relax. Discr ', DiscLevel{dl}), 'interpreter', 'latex', 'Fontsize', fs);
        xlabel( tl, '$ t \; [s] $', 'interpreter', 'latex', 'Fontsize', fs);
        ylabel( tl, '$x/L_x$', 'Interpreter', 'latex', 'FontSize', fs ) ;
        
        % colormap( flipud(gray(256)) ) ;
    
        % clearvars alpha_rho1 alpha_rho2 mom1 vel1 E alpha_rho_e1 alpha_rho_e2 pres tCoord xCoord
        
        % savefig(fig, fullfile(FigFolder, strcat(RelMod{rm}, DiscLevel{dl} ) ), '-v7.3' );
        % saveas(fig, fullfile(FigFolder, strcat(RelMod{rm}, DiscLevel{dl} ) ), 'epsc' );
        saveas(fig, fullfile(FigFolder, strcat(RelMod{rm}, DiscLevel{dl} ) ), 'png' );

        close
    end
end