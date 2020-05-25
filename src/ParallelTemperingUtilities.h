
/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef ParallelTemperingUtilities_H
#define ParallelTemperingUtilities_H

#include "GOMC_Config.h"    //For version number
#if GOMC_LIB_MPI
#include <mpi.h>
#endif
#include "XYZArray.h"
#include "ParallelTemperingPreprocessor.h"
#include "System.h"

class ParallelTemperingUtilities
{
public:

#if GOMC_LIB_MPI
explicit ParallelTemperingUtilities(MultiSim const*const& multisim, System & sys, StaticVals const& statV, ulong parallelTempFreq);
bool evaluateExchangeCriteria(ulong step);
void exchangePositions(XYZArray & myPos, MultiSim const*const& multisim);
void exchangeCOMs(XYZArray & myCOMs, MultiSim const*const& multisim);

#endif

private:

MultiSim const*const& ms;
SystemPotential & sysPotRef;
ulong parallelTempFreq;
vector<double> global_betas;

#if BOX_TOTAL == 1    
    vector<double> global_energies;
#else
    vector<vector <double>> global_energies;
#endif

};

#endif /*ParallelTemperingUtilities_H*/
