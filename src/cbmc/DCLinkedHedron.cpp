#define _USE_MATH_DEFINES
#include <math.h>
#include "DCLinkedHedron.h"
#include "DCData.h"
#include "TrialMol.h"
#include "../MolSetup.h"
#include "../Forcefield.h"
#include "../PRNG.h"
#include "../../lib/NumLib.h"
#include <numeric>
#include <cassert>

namespace
{
   struct FindA1
   {
      FindA1(uint x) : x(x) {};
      bool operator()(const mol_setup::Bond& b) { return (b.a1 == x); }
      uint x;
   };

   struct FindDih
   {
      FindDih(uint x, uint y) : x(x), y(y) {}
      uint x, y;
      bool operator()(const mol_setup::Dihedral d)
      { return (d.a0 == x && d.a3 == y) || (d.a0 == y && d.a3 == x); }
   };

}

namespace cbmc
{
      DCLinkedHedron::DCLinkedHedron
      (DCData* data,const mol_setup::MolKind& kind, uint focus, uint prev)
         : data(data), hed(data, kind, focus, prev)
      {
         using namespace mol_setup;
         using namespace std;
         vector<Bond> onPrev = AtomBonds(kind, hed.Prev());
         onPrev.erase(remove_if(onPrev.begin(), onPrev.end(), FindA1(hed.Focus())), onPrev.end());
         nPrevBonds = onPrev.size();

         for(uint i = 0; i < nPrevBonds; ++i)
	 {
            prevBonded[i] = onPrev[i].a1;
	    prevBondedLength[i] = data->ff.bonds.Length(onPrev[i].kind);
         }

	 vector<Bond> onFocus = AtomBonds(kind, hed.Focus());
         for(uint i = 0; i < onFocus.size(); ++i) {
            if (onFocus[i].a1 == prev) {
               focusPrevLength = data->ff.bonds.Length(onFocus[i].kind);
               break;
            }
         }

         vector<Dihedral> dihs = DihsOnBond(kind, hed.Focus(), hed.Prev());
         for(uint i = 0; i < hed.NumBond(); ++i)
	 {
            for(uint j = 0; j < nPrevBonds; ++j)
	    {
	       vector<Dihedral>::const_iterator match =
		 find_if(dihs.begin(), dihs.end(), FindDih(hed.Bonded(i),
							   prevBonded[j]));
               assert(match != dihs.end());
               dihKinds[i][j] = match->kind;
            }
         }
      }


   void DCLinkedHedron::PrepareNew(TrialMol& newMol, uint molIndex)
   {
     hed.PrepareNew(newMol, molIndex);
   }

   void DCLinkedHedron::PrepareOld(TrialMol& oldMol, uint molIndex)
   {
     hed.PrepareOld(oldMol, molIndex);
     for(uint i = 0; i < nPrevBonds; ++i)
     {
       prevBondedLengthOld[i] = sqrt(oldMol.OldDistSq(hed.Prev(), prevBonded[i]));
     }
       focusPrevLengthOld = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Prev()));
   }


   void DCLinkedHedron::BuildNew(TrialMol& newMol, uint molIndex)
   {
      PRNG& prng = data->prng;
      const CalculateEnergy& calc = data->calc;
      const Forcefield& ff = data->ff;
      uint nLJTrials = data->nLJTrialsNth;
      uint nDihTrials = data->nDihTrials;
      double* torsion = data->angles;
      double* torWeights = data->angleWeights;
      double* torEnergy = data->angleEnergy;
      double* ljWeights = data->ljWeights;
      double* bondedEn = data->bonded;
      double* inter = data->inter;
      double* nonbonded = data->nonbonded;
      double* nonbonded_1_4 = data->nonbonded_1_4;
      double* real = data->real;
      double* self = data->self;
      double* correction = data->correction;
      double* oneFour = data->oneFour;

      std::fill_n(ljWeights, nLJTrials, 0.0);
      std::fill_n(bondedEn, nLJTrials, 0.0);
      std::fill_n(oneFour, nLJTrials, 0.0);

      //get info about existing geometry
      newMol.SetBasis(hed.Focus(), hed.Prev());
      const XYZ center = newMol.AtomPosition(hed.Focus());
      XYZArray* positions = data->multiPositions;
      double prevPhi[MAX_BONDS];
      for (uint i = 0; i < hed.NumBond(); ++i)
      {
         //get position and shift to origin
         positions[i].Set(0, newMol.RawRectCoords(hed.BondLength(i),
						  hed.Theta(i), hed.Phi(i)));
      }
      for (uint i = 0; i < nPrevBonds; ++i)
      {
         double th;
         //not using theta, so this is a wasted cos and sqrt
         newMol.OldThetaAndPhi(prevBonded[i], hed.Prev(), th, prevPhi[i]);
      }
      XYZ rotationAxis = newMol.AtomPosition(hed.Focus()) -
	newMol.AtomPosition(hed.Prev());
      rotationAxis = data->axes.MinImage(rotationAxis, newMol.GetBox());
      rotationAxis *= (1 / rotationAxis.Length());
      RotationMatrix cross = RotationMatrix::CrossProduct(rotationAxis);
      RotationMatrix tensor = RotationMatrix::TensorProduct(rotationAxis);

      //counting backward to preserve prototype
      for (uint lj = nLJTrials; lj-- > 0;)
      {
	 ChooseTorsionNew(newMol, molIndex, prevPhi);
         ljWeights[lj] = std::accumulate(torWeights,
					 torWeights + nDihTrials, 0.0);
         uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
         bondedEn[lj] = torEnergy[winner];
	 oneFour[lj] = nonbonded_1_4[winner];
         //convert chosen torsion to 3D positions
         RotationMatrix spin = RotationMatrix::FromAxisAngle(-torsion[winner],
							     cross, tensor);
         for (uint b = 0; b < hed.NumBond(); ++b)
	 {
               //find positions
               positions[b].Set(lj, spin.Apply(positions[b][0]));
               positions[b].Add(lj, center);
         }
      }

      for (uint b = 0; b < hed.NumBond(); ++b)
      {
         data->axes.WrapPBC(positions[b], newMol.GetBox());
      }

      double stepWeight = EvalLJ(newMol, molIndex);
      uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
      for(uint b = 0; b < hed.NumBond(); ++b)
      {
         newMol.AddAtom(hed.Bonded(b), positions[b][winner]);
      }
      newMol.AddEnergy(Energy(bondedEn[winner] + hed.GetEnergy(),
			      nonbonded[winner] + hed.GetNonBondedEn() +
			      oneFour[winner], inter[winner], real[winner],
			      0.0, self[winner], correction[winner]));
      newMol.MultWeight(hed.GetWeight());
      newMol.MultWeight(stepWeight);
   }

   void DCLinkedHedron::BuildOld(TrialMol& oldMol, uint molIndex)
   {
      PRNG& prng = data->prng;
      const CalculateEnergy& calc = data->calc;
      const Forcefield& ff = data->ff;
      uint nLJTrials = data->nLJTrialsNth;
      uint nDihTrials = data->nDihTrials;
      double* torsion = data->angles;
      double* torWeights = data->angleWeights;
      double* torEnergy = data->angleEnergy;
      double* ljWeights = data->ljWeights;
      double* bondedEn = data->bonded;
      double* inter = data->inter;
      double* nonbonded = data->nonbonded;
      double* nonbonded_1_4 = data->nonbonded_1_4;
      double* real = data->real;
      double* self = data->self;
      double* correction = data->correction;
      double* oneFour = data->oneFour;

      std::fill_n(ljWeights, nLJTrials, 0.0);
      std::fill_n(bondedEn, nLJTrials, 0.0);
      std::fill_n(oneFour, nLJTrials, 0.0);
      std::fill_n(oneFour, nLJTrials, 0.0);

      //get info about existing geometry
      oldMol.SetBasis(hed.Focus(), hed.Prev());
      //Calculate OldMol Bond Energy &
      //Calculate phi weight for nTrials using actual theta of OldMol
      hed.ConstrainedAnglesOld(data->nAngleTrials - 1, oldMol, molIndex);
      const XYZ center = oldMol.AtomPosition(hed.Focus());
      XYZArray* positions = data->multiPositions;
      double prevPhi[MAX_BONDS];
      for (uint i = 0; i < hed.NumBond(); ++i)
      {
         //get position and shift to origin
         positions[i].Set(0, oldMol.AtomPosition(hed.Bonded(i)));
         data->axes.UnwrapPBC(positions[i], 0, 1, oldMol.GetBox(), center);
         positions[i].Add(0, -center);
      }
      for (uint i = 0; i < nPrevBonds; ++i)
      {
         double t;
         //not using theta, so this is a wasted cos and sqrt
         oldMol.OldThetaAndPhi(prevBonded[i], hed.Prev(), t, prevPhi[i]);
      }
      XYZ rotationAxis = oldMol.AtomPosition(hed.Focus()) - 
      oldMol.AtomPosition(hed.Prev());
      rotationAxis = data->axes.MinImage(rotationAxis, oldMol.GetBox());
      rotationAxis *= (1 / rotationAxis.Length());
      RotationMatrix cross = RotationMatrix::CrossProduct(rotationAxis);
      RotationMatrix tensor = RotationMatrix::TensorProduct(rotationAxis);

      //counting backward to preserve prototype
      for (uint lj = nLJTrials; lj-- > 1;)
      {
	 ChooseTorsionOld(oldMol, molIndex, prevPhi);
         ljWeights[lj] = std::accumulate(torWeights, torWeights + nDihTrials,
					 0.0);
         uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
         bondedEn[lj] = torEnergy[winner];
	 oneFour[lj] = nonbonded_1_4[winner];
         //convert chosen torsion to 3D positions
         RotationMatrix spin =
	   RotationMatrix::FromAxisAngle(-torsion[winner], cross, tensor);
         for (uint b = 0; b < hed.NumBond(); ++b)
	 {
               //find positions
               positions[b].Set(lj, spin.Apply(positions[b][0]));
               positions[b].Add(lj, center);
         }
      }
      ljWeights[0] = 0.0;
      for (uint tor = 0; tor < nDihTrials; ++tor)
      {
         torsion[tor] = (tor == 0) ? 0.0 : data->prng.rand(M_PI * 2);
         torEnergy[tor] = 0.0;
	 nonbonded_1_4[tor] = 0.0;
         for (uint b = 0; b < hed.NumBond(); ++b)
	 {
	    double theta1 =  hed.Theta(b);
            double trialPhi = hed.Phi(b) + torsion[tor];
            for (uint p = 0; p < nPrevBonds; ++p)
	    {
	       double theta0 = oldMol.GetTheta(prevBonded[p], hed.Prev(),
					    hed.Focus());
	       double distSq =
		 oldMol.DihedDist(prevBondedLengthOld[p], focusPrevLengthOld,
				  hed.BondLengthOld(b), theta0, theta1,
				  trialPhi - prevPhi[p]);
	       nonbonded_1_4[tor] +=
		 data->calc.IntraEnergy_1_4(distSq, prevBonded[p],
					    hed.Bonded(b), molIndex);
	       if(isnan(nonbonded_1_4[tor]))
		 nonbonded_1_4[tor] = num::BIGNUM;

               torEnergy[tor] += ff.dihedrals.Calc(dihKinds[b][p],
						   trialPhi - prevPhi[p]);
            }
         }
         ljWeights[0] += exp(-ff.beta * (torEnergy[tor] + nonbonded_1_4[tor]));
      }
      bondedEn[0] = torEnergy[0];
      oneFour[0] = nonbonded_1_4[0];

      for (uint b = 0; b < hed.NumBond(); ++b)
      {
         positions[b].Add(0, center);
         data->axes.WrapPBC(positions[b], oldMol.GetBox());
      }
      double stepWeight = EvalLJ(oldMol, molIndex);
      for(uint b = 0; b < hed.NumBond(); ++b)
      {
         oldMol.ConfirmOldAtom(hed.Bonded(b));
      }
      oldMol.AddEnergy(Energy(bondedEn[0] + hed.GetEnergy() +
			      hed.GetOldBondEn(), nonbonded[0] +
			      hed.GetNonBondedEn() + oneFour[0],
			      inter[0], real[0], 0.0, self[0],
			      correction[0]));
 
      oldMol.MultWeight(hed.GetWeight());
      oldMol.MultWeight(stepWeight);
   }

   double DCLinkedHedron::EvalLJ(TrialMol& mol, uint molIndex)
   {
      uint nLJTrials = data->nLJTrialsNth;
      double* inter = data->inter;
      double* nonbonded = data->nonbonded;
      //double* nonbonded_1_4 = data->nonbonded_1_4;
      double* real = data->real;
      double* self = data->self;
      double* correction = data->correction;
      XYZArray* positions = data->multiPositions;

      std::fill_n(data->inter, nLJTrials, 0.0);
      std::fill_n(data->nonbonded, nLJTrials, 0.0);
      //std::fill_n(nonbonded_1_4, nLJTrials, 0.0);
      std::fill_n(self, nLJTrials, 0.0);
      std::fill_n(real, nLJTrials, 0.0);
      std::fill_n(correction, nLJTrials, 0.0);

       for (uint b = 0; b < hed.NumBond(); ++b)
      {
	data->calc.ParticleInter(inter, real, positions[b], hed.Bonded(b),
                                  molIndex, mol.GetBox(), nLJTrials);
	data->calc.ParticleNonbonded(nonbonded, mol, positions[b],
				     hed.Bonded(b), mol.GetBox(),
				     nLJTrials);
<<<<<<< HEAD
	data->calcEwald.SwapSelf(self, molIndex, hed.Bonded(b), mol.GetBox(),
				 nLJTrials);
	data->calcEwald.SwapCorrection(correction, mol, positions, b, 
				       hed.bonded, mol.GetBox(), nLJTrials,
				       hed.Prev(), false);
=======
>>>>>>> development
	//data->calc.ParticleNonbonded_1_4(nonbonded_1_4, mol, positions[b], 
	//				 hed.Bonded(b), mol.GetBox(),
	//				 nLJTrials);

      }
      double stepWeight = 0;
      for (uint lj = 0; lj < nLJTrials; ++lj)
      {
	 data->ljWeights[lj] *= exp(-data->ff.beta * 
				    (inter[lj] + nonbonded[lj] + real[lj] +
				     self[lj] + correction[lj]));
         stepWeight += data->ljWeights[lj];
      }
      return stepWeight;
   }

   void DCLinkedHedron::ChooseTorsionNew(TrialMol& newMol, uint molIndex,
					 double prevPhi[])
   {
      double* torsion = data->angles;
      double* torEnergy = data->angleEnergy;
      double* torWeights = data->angleWeights;
      double* nonbonded_1_4 = data->nonbonded_1_4;
      uint nDihTrials = data->nDihTrials;
      const Forcefield& ff = data->ff;

      std::fill_n(torsion, data->nDihTrials, 0.0);
      std::fill_n(torWeights, data->nDihTrials, 0.0);
      std::fill_n(torEnergy, data->nDihTrials, 0.0);
      std::fill_n(nonbonded_1_4, data->nDihTrials, 0.0);

      //select torsion based on all dihedral angles
      for (uint tor = 0; tor < nDihTrials; ++tor)
      {
         torsion[tor] = data->prng.rand(M_PI * 2);
         torEnergy[tor] = 0.0;
	 nonbonded_1_4[tor] = 0.0;
         for (uint b = 0; b < hed.NumBond(); ++b)
	 {
	    double theta1 =  hed.Theta(b);
            double trialPhi = hed.Phi(b) + torsion[tor];
            for (uint p = 0; p < nPrevBonds; ++p)
	    {
	       double theta0 = newMol.GetTheta(prevBonded[p], hed.Prev(),
					    hed.Focus());
	       double distSq =
		 newMol.DihedDist(prevBondedLength[p], focusPrevLength,
			       hed.BondLength(b), theta0, theta1,
			       trialPhi - prevPhi[p]);
	       nonbonded_1_4[tor] +=
		 data->calc.IntraEnergy_1_4(distSq, prevBonded[p],
					    hed.Bonded(b), molIndex);
	       if(isnan(nonbonded_1_4[tor]))
		 nonbonded_1_4[tor] = num::BIGNUM;

               torEnergy[tor] += ff.dihedrals.Calc(dihKinds[b][p],
						   trialPhi - prevPhi[p]);
            }
	 }        
	 torWeights[tor] = exp(-ff.beta *(torEnergy[tor] + nonbonded_1_4[tor]));   
      }
   }

   void DCLinkedHedron::ChooseTorsionOld(TrialMol& oldMol, uint molIndex,
					 double prevPhi[])
   {
      double* torsion = data->angles;
      double* torEnergy = data->angleEnergy;
      double* torWeights = data->angleWeights;
      double* nonbonded_1_4 = data->nonbonded_1_4;
      uint nDihTrials = data->nDihTrials;
      const Forcefield& ff = data->ff;

      std::fill_n(torsion, data->nDihTrials, 0.0);
      std::fill_n(torWeights, data->nDihTrials, 0.0);
      std::fill_n(torEnergy, data->nDihTrials, 0.0);
      std::fill_n(nonbonded_1_4, data->nDihTrials, 0.0);

      //select torsion based on all dihedral angles
      for (uint tor = 0; tor < nDihTrials; ++tor)
      {
         torsion[tor] = data->prng.rand(M_PI * 2);
         torEnergy[tor] = 0.0;
	 nonbonded_1_4[tor] = 0.0;
         for (uint b = 0; b < hed.NumBond(); ++b)
	 {
	    double theta1 =  hed.Theta(b);
            double trialPhi = hed.Phi(b) + torsion[tor];
            for (uint p = 0; p < nPrevBonds; ++p)
	    {
	       double theta0 = oldMol.GetTheta(prevBonded[p], hed.Prev(),
					    hed.Focus());
	       double distSq =
		 oldMol.DihedDist(prevBondedLengthOld[p], focusPrevLengthOld,
			       hed.BondLengthOld(b), theta0, theta1,
			       trialPhi - prevPhi[p]);
	       nonbonded_1_4[tor] +=
		 data->calc.IntraEnergy_1_4(distSq, prevBonded[p],
					    hed.Bonded(b), molIndex);
	       if(isnan(nonbonded_1_4[tor]))
		 nonbonded_1_4[tor] = num::BIGNUM;

               torEnergy[tor] += ff.dihedrals.Calc(dihKinds[b][p],
						   trialPhi - prevPhi[p]);
            }
	 }        
	 torWeights[tor] = exp(-ff.beta *(torEnergy[tor] + nonbonded_1_4[tor]));
      }
   }
   
}
