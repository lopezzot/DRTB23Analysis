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

// Return true for e+ using PreShower and MuonCounter only
bool IsPositronPsMu(const double& PreShower, const double& MuonTrk)
{
  const int PScut = 1500;
  const int MuonTrkcut = 200;
  if (PreShower > PScut && MuonTrk < MuonTrkcut)
    return true;
  else
    return false;
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

// Return true is total S SiPM energy is above cut
bool IsSiPMSabovecut(const float svec[160], const double& cut)
{
  return std::accumulate(svec, svec + 160, 0.) > cut ? true : false;
}

// Return true if event is in DWC radius
bool IsDWCradius(double pos[2], const double& radiuscut, const int& DWCIdx)
{
  double xoffset = 0.;
  double yoffset = 0.;
  if (DWCIdx == 1) {  // DWC 1
    xoffset = 1.17;
    yoffset = -4.3;
  }
  else if (DWCIdx == 2) {  // DWC 2
    xoffset = 4.3;
    yoffset = 0.0;
  }
  else {
    std::cout << "Wrong DWC index (1 or 2), going to std::abort()." << std::endl;
    std::abort();
  }
  pos[0] = pos[0] + xoffset;
  pos[1] = pos[1] + yoffset;
  double radius = std::sqrt(pow(pos[0], 2.) + pow(pos[1], 2.));

#ifdef DEBUG
  std::cout << "DWC" << DWCIdx << " radius " << radius << " mm." << std::endl;
#endif

  return (radius < radiuscut) ? true : false;
}

}  // namespace utils

//**************************************************