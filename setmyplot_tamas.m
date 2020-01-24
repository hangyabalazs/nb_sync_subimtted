function  setmyplot_tamas(ax,varargin)
%SETMYPLOT_TAMAS   Set axis properties.
%   SETMYPLOT_TAMAS(AX) sets axis properties.
%
%   See also SETMYPLOT and SETMYFIGURE.


% Input arguments
if nargin == 0
    ax = gca;
end

% Set axis properties
set(ax,'TickDir','out','box','off');
set(ax,'FontSize',24,'LineWidth',1);
set(ax,'TickLength',[0.05 0])
set(ax,'FontName','Arial');

% Title and axis label font size
th = findobj(allchild(gca),'Type','text','VerticalAlignment','bottom','Rotation',0);
axlab = setdiff(findobj(allchild(gca),'Type','text'),th);
if ~isempty(th)
    set(th,'FontSize',24,'FontName','Arial')
end
if ~isempty(axlab),
    set(axlab,'FontSize',24,'FontName','Arial')
end
if nargin > 1
    set(varargin{:});
end