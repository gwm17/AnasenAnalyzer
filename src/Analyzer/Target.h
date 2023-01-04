#ifndef TARGET_H
#define TARGET_H

#include <string>
#include <vector>
#include <cmath>
#include "catima/gwm_integrators.h"
#include "MassLookup.h"

namespace AnasenAnalyzer {

	class Target {
	
	public:
	 	Target() = default;
	 	Target(const std::vector<uint32_t>& z, const std::vector<uint32_t>& a, const std::vector<int>& stoich, double density);
	 	~Target();

	 	double GetEnergyLoss(int zp, int ap, double startEnergy, double pathLength);
	 	double GetReverseEnergyLoss(int zp, int ap, double finalEnergy, double pathLength);
		double GetPathLength(int zp, int ap, double startEnergy, double finalEnergy); //Returns pathlength for a particle w/ startE to reach finalE (cm)
		double GetAngularStraggling(int zp, int ap, double energy, double pathLength); //Returns planar angular straggling in radians for a particle with energy and pathLength
	 	inline double GetDensity() { return m_material.density(); } //g/cm^3
	
	private:
		catima::Material m_material;

		static constexpr double s_epsilon = 1.0e-6;
	};

}

#endif 