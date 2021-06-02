/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FFImpropers.h"  //Parent class
#include "FFSetup.h" //For initialization data
#include <algorithm> //for vector copying


void FFImpropers::Init(ff_setup::Improper const& imp)
{
  uint size = dih.getTerms(), numSubDiv = dih.getnamecnt(), count = 0;
  Kchi = new double[size];
  n = new uint[size];
  delta = new double[size];
  subdiv.Init(numSubDiv);
  for (uint s = 0; s < numSubDiv; s++) {
    std::string div = dih.getname(s);
    uint cnt = dih.append(div, Kchi, delta, n, count);
    subdiv.Set(s, count, cnt);
    count += cnt;
  }
}

FFImpropers::~FFImpropers()
{
  delete[] Komega;
  delete[] omega0;
  delete[] n;
}
