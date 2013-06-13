
function plotFigure(varargin)

figColor = 'b';
figMarker = 'p';
figLineType = '-';
figLineWidth = 1;

xValues = varargin{1,1};
yValues = varargin{1,2};
figPage = varargin{1,3};
figType = varargin{1,4};

figure(figPage);hold all;grid on;

switch figType
    
    case 'plot'
        
        plot(xValues,yValues,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
        
    case 'semilogy'
        
        semilogy(xValues,yValues,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);
        
    case 'cdfplot'
        
        cdfplot(xValues,'Color',figColor,'LineWidth',figLineWidth,...
            'LineStyle',figLineType,'MarkerFaceColor',figColor,'Marker',figMarker);

end
        