/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_IMPROPERS_H
#define FF_IMPROPERS_H

#include <cmath>       //cos, pow

#include "FFConst.h"    //GetRot
#include "SubdividedArray.h" //Subdivisions in master array for dih sets.

#include "NumLib.h"    //Sq
#include "BasicTypes.h" //For "uint"

namespace ff_setup
{
class Improper;
}

//FFDihKind
//Stores parameters for Impropers kinds
//dihedrals may require several sets of parameters to calculate,
//these are stored in consecutive array elements, each kind's
//parameter sets are delinated by subdiv
//Dihedrals may be calculated by several different styles
class FFImpropers
{
public:
  //calculate the energy of dih kind at angle phi
  double Calc(const uint kind, const double phi) const;
  //Initialize with data from parameter files
  void Init(ff_setup::Improper const& imp);

  FFImpropers(void) : Kchi(NULL), delta(NULL), n(NULL) {}
  ~FFImpropers(void);

private:

  //dih kind params
  SubdividedArray subdiv;
  double * Komega, * omega0;
  uint *n;
};

inline double FFDihedrals::Calc(const uint kind, const double phi) const
{
  double sum = 0.0;
  for(uint i = subdiv.Begin(kind); i != subdiv.End(kind); ++i) {
    sum += Kchi[i] * (1 + cos(n[i] * phi - delta[i]));
  }
  return sum;
}

#endif /*FF_DIHEDRALS_H*/
