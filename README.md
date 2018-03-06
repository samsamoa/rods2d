Simulation of Self-Propelled Rods in 2D
=======================================

This code was used for the academic paper "Spontaneous Segregation of Self-Propelled Particles with Different Motilities", published in [Soft Matter (Jan 2012)](http://pubs.rsc.org/en/Content/ArticleLanding/2012/SM/c2sm06960a) and available open-access at [arXiv](https://arxiv.org/abs/1110.2479).

As is typical for such projects, this code was used only by myself, for only this paper.  I have made no attempt to clean it up or add additional documentation (other than this ReadMe file).

The key components are:
- **RodSimulation**: Keeps track of the components and parameters of the simulation 
- **Rod**: A propelled rod (the primary element of the simulation)
- **InteractionPair**: For efficiency, we only apply interaction forces to a pair of rods when they are close enough to interact.  This class keeps track of a pair that is close enough to potentially interact.

In addition, there is some code used for data analysis (clumps.cpp, densityFluctuations.cpp, and the scripts/ directory).

Pictures and videos of the results can be found in the paper itself, and on the sidebar of the [arXiv](https://arxiv.org/abs/1110.2479) page.
