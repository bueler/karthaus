Here is a project suggestion for 2014:

HOW FAST DOES A REGIONAL MODEL DIVERGE FROM A SURROUNDING WHOLE ICE SHEET MODEL?
--------------------------------------------------------------------------------

In the lectures we already have a (minimal) isothermal SIA model for the whole Antarctic ice sheet.  Run this first; see `mfiles/ant/ant.m` and `mfiles/ant/buildant.m`.  It is at 50 km resolution.  Also build the data for 25 km resolution, and run that.  Note the performance difference.

The achievable resolution of such a model is limited by the large size of the ice sheet.  The (floating point) computational time in the stress balance and the memory access time, in particular, limit the ability to do large ice sheet computations.  Also, stability properties of the mass continuity scheme determine how many time steps must be computed.  If the full three-dimensional temperature or velocity fields are stored, the achievable resolution can be limited by the raw amount of memory or even hard disk space available.

For reasons like these, one might build a _regional_ model with higher resolution but having smaller spatial extent.  Also, data may be available for the region that is not available for the whole ice sheet; perhaps an individual outlet glacier or ice stream/shelf has been more densely observed.  Because regional models do not cover the whole ice sheet, they need boundary conditions connecting them to the rest of the ice sheet.

I propose that you set up a regional model for a part of Antarctica, starting from my whole ice sheet setup.  You can compare the volume and flow speed within a region _R_ of a whole ice sheet model to the same properties from the regional model of _R_ alone.  That is, there are two numerical models that cover _R_, one being the whole ice sheet model and the other just covering _R_, the regional model.  The boundary conditions for the second (regional) model are provided initially and then held fixed in time.  In the simplest comparison, which you should do, the grid resolution is the same in both models, so that the only difference is that nontrivial thickness boundary conditions apply in the regional model.  You will run the two models for the same period and using the same starting values (i.e.~of ice thickness) and compare the results.

You choose the region to model.  You choose how the boundary conditions come from the whole ice sheet model.  You will need to think about how to cut out the region, and how to modify the main model code (i.e.~`mfiles/siageneral.m`) to use fixed thickness boundary conditions.

The issue is to understand the regional "mode" of modeling, and its advantages and disadvantages, so the flow model should be the same isothermal SIA choice that is already implemented.  Numerically approximating new (continuum) model equations is not the issue!

References/examples for this general topic:

  * I. Joughin, B. Smith, and B. Medley (2014). _Marine ice sheet collapse potentially under way for the Thwaites Glacier Basin, West Antarctica_, Science 344 (6185), 735--738, doi:10.1126/science.1249055.
  * M. Koutnick and E. Waddington (2012).  _Well-posed boundary conditions for limited-domain models of transient ice flow near an ice divide_, J. Glaciol. 58 (211), 1008--1020, doi: 10.3189/2012JoG11J212.
  * M. Mengel and A. Levermann (2014), _Ice plug prevents irreversible discharge from East Antarctica_, Nature Clim. Change, doi: 10.1038/nclimate2226.

