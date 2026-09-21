set(groot,'defaultFigureCreateFcn',@(fig, ~)addFigButtons(fig));
function addFigButtons(fig)
  if ~matlab.ui.internal.isUIFigure(fig)
     addToolbarExplorationButtons(fig)
  end
end