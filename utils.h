//**************************************************
// \file utils.h
// \brief: Utils for Analysis of 20223 TB
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
// 	    @lopezzot
// \start date: 6 July 2023
//**************************************************

#include "PhysicsEvent.h"

#include <cstdlib>

// #define DEBUG

namespace utils
{

// Return true for preshower values aroung 8000
bool IsPSUltraHigh(const double& PreShower)
{
  // in the 2023 test-beam data there are few events with
  // preshower values aroung 8000.
  // It seems a software artifact, so I remove them
  if (PreShower > 7000.)
    return true;
  else
    return false;
}

// Return true for e+ using PreShower and MuonCounter only
bool IsPositronPsMu(const double& PreShower, const double& MuonTrk)
{
  const int PScut = 700;
  const int MuonTrkcut = 200;
  if (PreShower > PScut && MuonTrk < MuonTrkcut)
    return true;
  else
    return false;
}

// Return true is PreShower is above cut value
bool IsPsAboveCut(const double& PreShower, const double& cut)
{
  if (PreShower > cut)
    return true;
  else
    return false;
}

// Return true for pi+ using PreShower and MuonCounter only
bool IsPionPsMu(const double& PreShower, const double& MuonTrk)
{
  const int PScut = 700;
  const int PSpedestalcut = 500;
  const int MuonTrkcut = 200;
  if (PreShower > PSpedestalcut && PreShower < PScut && MuonTrk < MuonTrkcut)
    return true;
  else
    return false;
}

// Returns true if pion has interacted in the central tower (dummy right now)
bool HasPionInteracted(const double& Senergy, const double& Pioncutenergy)
{
  return (Senergy > Pioncutenergy) ? true : false;
}

// Return SiPM S position from index
double* ScinSiPMmap(const int& index)
{
  static double SSiPMpos[2];
  int row = index / 16;
  int column = (index - 16 * row);
  SSiPMpos[0] = -1.5 + 2 * (column - 7);
  SSiPMpos[1] = 2. * sqrt(3.) * (4 - row) + sqrt(3) / 2.;
  return SSiPMpos;
}

// Return SiPM C position form index
double* CherSiPMmap(const int& index)
{
  static double CSiPMpos[2];
  int row = index / 16;
  int column = (index - 16 * row);
  CSiPMpos[0] = -0.5 + 2.0 * (column - 7);
  CSiPMpos[1] = 1.5 * sqrt(3) + 2. * sqrt(3) * (4 - row);

  return CSiPMpos;
}

// Return shower barycenter from vector of S SiPMS
double* GetScinbar(const float svec[160])
{
  static double Sbar[2];
  double x = 0;
  for (unsigned int index = 0; index < 160; index++) {
    x += svec[index] * ScinSiPMmap(index)[0];
  }
  x = x / std::accumulate(svec, svec + 160, 0.);

  double y = 0;
  for (unsigned int index = 0; index < 160; index++) {
    y += svec[index] * ScinSiPMmap(index)[1];
  }
  y = y / std::accumulate(svec, svec + 160, 0.);
  Sbar[0] = x;
  Sbar[1] = y;

  return Sbar;
}

// Return shower barycenter from vector of C SiPMS
double* GetCherbar(const float cvec[160])
{
  static double Cbar[2];
  double x = 0;
  for (unsigned int index = 0; index < 160; index++) {
    x += cvec[index] * CherSiPMmap(index)[0];
  }
  x = x / std::accumulate(cvec, cvec + 160, 0.);

  double y = 0;
  for (unsigned int index = 0; index < 160; index++) {
    y += cvec[index] * CherSiPMmap(index)[1];
  }
  y = y / std::accumulate(cvec, cvec + 160, 0.);
  Cbar[0] = x;
  Cbar[1] = y;

  return Cbar;
}

// Return distance between two (x,y) points (e.g. fiber distance to baricenter)
double Getdist(double pos[2], double bar[2])
{
  double radius = std::sqrt(pow(pos[0] - bar[0], 2.) + pow(pos[1] - bar[1], 2.));

  return radius;
}

// Return true is total S SiPM energy is above cut
bool IsSiPMSabovecut(const float svec[160], const double& cut)
{
  return std::accumulate(svec, svec + 160, 0.) > cut ? true : false;
}

// Return true if event is in DWC radius
bool IsDWCradius(double pos[2], const double& radiuscut, const std::array<double, 2>& offset)
{
  double newpos[2]{};
  newpos[0] = pos[0] + offset[0];
  newpos[1] = pos[1] + offset[1];
  double radius = std::sqrt(pow(newpos[0], 2.) + pow(newpos[1], 2.));

#ifdef DEBUG
  std::cout << "DWC"
            << " radius " << radius << " mm." << std::endl;
#endif

  return (radius < radiuscut) ? true : false;
}

// Return true if event is in DWC radius, no offset expected
bool IsDWCradius(double pos[2], const double& radiuscut)
{
  double radius = std::sqrt(pow(pos[0], 2.) + pow(pos[1], 2.));

#ifdef DEBUG
  std::cout << "DWC"
            << " radius " << radius << " mm." << std::endl;
#endif

  return (radius < radiuscut) ? true : false;
}

}  // namespace utils

//**************************************************
