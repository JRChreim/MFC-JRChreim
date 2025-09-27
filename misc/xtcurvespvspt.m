% Example of using binary_reader_wrapper
% Need to change "format" to 2 when post-processing the example cases

clear; clc; close all;

% toggle location
loc = 'local';
switch loc
    case 'carpenter'
        mfcPath = '/p/global/jrchreim/simulations/PhaseChange/1D/BubbleDynamics/StrongCollapse/6Eqn/';
        RelMod = {'pTFinalTest', 'pFinalTest'} ;
        % RelMod = {'pTFinalTest'} ;
        DiscLevel = {'N1280E3'} ;
        % DiscLevel = {'N160E3'} ;
        compliment = 'Cartesian/BC-6/C000E-00';
        FigFolder = '/p/global/jrchreim/Figures/';
    case 'local'
        mfcPath = '/disk/simulations/PhaseChange/ShockTube/1D/StrongCollapse/6Eqn/';
        RelMod = {'p'} ;
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
   
        % normalizing variables
        Lx = xCoord(end,1) - xCoord(1,1) ;        

        tOtend = tCoord ./ tCoord(end) ;
        xOL = xCoord ./ Lx ;

        % filtering out information based on the presence of vapor
        % (focusing on the vapor region)

        % SETTING UP minimum and maximum values for the colorbars
        pres(alpha2 < 1 - 1E-8) = NaN ;
        T1(alpha2 < 1 - 1E-8) = NaN ;
        T2(alpha2 < 1 - 1E-8) = NaN ;
        alpha2(alpha2 < 1 - 1E-8) = NaN ;
        if ( dl == 1 && rm == 1)
            maxP  = max( max( pres ) ) ; minP = min( min( pres ) ) ;
            maxTv = max( max( T2 ) )  ; minTv = min( min( T2 ) ) ;
            maxTl = max( max( T1 ) ) ; minTl = min( min( T1 ) ) ;
            maxM2 = max( max( mom1 ) ) ; minM2 = min( min( mom1 ) ) ;
        end

        % setting up figure
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        tl = tiledlayout( 2, 2, 'TileSpacing', 'compact' ); fs = 30 ;

        % PLOTTING VARIABLES
        % pressure
        % levP = floor(minP):ceil(maxP) ;
        levP = 8;

        xC{1} = (xOL >= 0.0 & xOL <= 7.0E-3) ;
        XM = reshape(tOtend(xC{1}), [], size(xOL, 2)) ; YM = reshape(xOL( xC{1} ), [], size(xOL, 2)) ;
        ZM = reshape(pres(xC{1}), [], size(xOL,2) ) ;

        nexttile(1) ; contourf(XM, YM, ZM, levP) ; colorbar ;
        title( '$ p \; [Pa] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.2f'); ytickformat('%.2f');
        ax = gca; ax.FontSize = fs;

        text(0.5, 0.30, ['p_{max}: ' num2str( max( max( pres ) ) ) ], 'FontSize', fs);
        text(0.5, 0.25, ['p_{min}: ' num2str( min( min( pres ) ) ) ], 'FontSize', fs);
         
        % Vapor Temperature
        % levTv = floor(minTv):10:ceil(maxTv) ;
        levTv = 8 ;
        xC{2} = (xOL >= 0.0 & xOL <= 7.0E-3) ;
        XM = reshape(tOtend(xC{2}), [], size(xOL, 2)) ; YM = reshape(xOL( xC{2} ), [], size(xOL, 2)) ;
        ZM = reshape(T2(xC{2}), [], size(xOL,2) ) ;

        nexttile(2) ; contourf(XM, YM, ZM, levTv) ; colorbar ;
        title( '$ T_{v} \; [K] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.2f');
        ax = gca; ax.FontSize = fs;

        text(0.5, 0.30, ['T_{max}: ' num2str( max( max( T2 ) ) ) ], 'FontSize', fs);
        text(0.5, 0.25, ['T_{min}: ' num2str( min( min( T2 ) ) ) ], 'FontSize', fs);

        % liquid temperature
        % levTl = floor(minTl):ceil(maxTl) ;
        levTl = 8 ;
        xC{3} = (xOL >= 0.0 & xOL <= 7.0E-3) ;
        XM = reshape(tOtend(xC{3}), [], size(xOL, 2)) ; YM = reshape(xOL( xC{3} ), [], size(xOL, 2)) ;
        ZM = reshape(T1(xC{3}), [], size(xOL,2) ) ;

        nexttile(3) ; contourf(XM, YM, ZM, levTl) ; colorbar ;
        title( '$ T_{l} \; [K] $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.2f');
        ax = gca; ax.FontSize = fs;

        text(0.5, 0.30, ['T_{max}: ' num2str( max( max( T1 ) ) ) ], 'FontSize', fs);
        text(0.5, 0.25, ['T_{min}: ' num2str( min( min( T1 ) ) ) ], 'FontSize', fs);

        % levTV4 = floor(minM2):ceil(maxM2) ;
        levTV4 = 10 ;
        xC{4} = (xOL >= 0.0 & xOL <= 7.0E-3) ;
        XM = reshape(tOtend(xC{4}), [], size(xOL, 2)) ; YM = reshape(xOL( xC{4} ), [], size(xOL, 2)) ;
        ZM = reshape(alpha2(xC{4}), [], size(xOL,2) ) ;

        nexttile(4) ; contourf(XM, YM, ZM, 'EdgeColor', 'none') ; colorbar ;
        title( '$ \alpha_{v} \;  $', 'interpreter', 'latex', 'Fontsize', fs);
        xtickformat('%.1f'); ytickformat('%.4f');
        ax = gca; ax.FontSize = fs;

        % text(0.5,0.325,['m_{2,max}: ' num2str( max( max( alpha_rho2 ) ) ) ], 'FontSize', fs);
        % text(0.5,0.335,['m_{2,min}: ' num2str( min( min( alpha_rho2 ) ) ) ], 'FontSize', fs);

        title( tl, strcat(RelMod{rm}, '-relax. Discr ', DiscLevel{dl}), 'interpreter', 'latex', 'Fontsize', fs);
        xlabel( tl, '$ t/t_{final} \;  $', 'interpreter', 'latex', 'Fontsize', fs);
        ylabel( tl, '$x/L_x$', 'Interpreter', 'latex', 'FontSize', fs ) ;

        clearvars alpha_rho1 alpha_rho2 mom1 vel1 E alpha_rho_e1 alpha_rho_e2 pres tCoord xCoord

        % savefig(fig, fullfile(FigFolder, strcat(RelMod{rm}, DiscLevel{dl} ) ), '-v7.3' );
        % saveas(fig, fullfile(FigFolder, strcat(RelMod{rm}, DiscLevel{dl} ) ), 'epsc' );
        % saveas(fig, fullfile(FigFolder, strcat(RelMod{rm}, DiscLevel{dl} ) ), 'png' );

        % close
    end
end