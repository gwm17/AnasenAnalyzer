#include "Target.h"
#include "catima/nucdata.h"
#include "Precision.h"

#include <iostream>

namespace AnasenAnalyzer {

	//z,a: istope list of material compound, stoich: compound stoichometry, density: material density in g/cm^3
	Target::Target(const std::vector<uint32_t>& z, const std::vector<uint32_t>& a, const std::vector<int>& stoich, double density)
	{
		MassLookup& masses = MassLookup::GetInstance();
		for(size_t i=0; i<z.size(); i++)
		{
			m_material.add_element(masses.FindMassU(z[i], a[i]), z[i], stoich[i]);
		}

		m_material.density(density); //g/cm^3
	}
	
	Target::~Target() {}
	
	/*Calculates energy loss for travelling all the way through the target*/
	//ZP, AP: projectile isotope, startEnergy: MeV, pathLength: m
	//return eloss in MeV
	double Target::GetEnergyLoss(int zp, int ap, double startEnergy, double pathLength)
	{
		if(Precision::IsFloatLessOrAlmostEqual(startEnergy, 0.0, s_epsilon) || std::isnan(pathLength))
			return 0.0;
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = startEnergy/proj.A;
		m_material.thickness_cm(pathLength * 100.0); //Takes in a path length and calculates density corrected thickness to thickness param
		return  startEnergy - (catima::energy_out(proj, m_material) * proj.A);
	}

	/*Calculates reverse energy loss for travelling all the way through the target*/
	//ZP, AP: projectile isotope, finalEnergy: MeV, pathLength: m
	//return eloss in MeV
	double Target::GetReverseEnergyLoss(int zp, int ap, double finalEnergy, double pathLength)
	{
		if(Precision::IsFloatLessOrAlmostEqual(finalEnergy, 0.0, s_epsilon) || std::isnan(pathLength))
			return 0.0;
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = finalEnergy/proj.A;
		m_material.thickness_cm(pathLength * 100.0);
		return catima::reverse_integrate_energyloss(proj, m_material);
	}

	//Get the path length (range) for a particle with incoming energy startEnergy and a outgoing energy finalEnergy 
	double Target::GetPathLength(int zp, int ap, double startEnergy, double finalEnergy)
	{
		double densityInv = 1.0/m_material.density();
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = startEnergy/proj.A;
		double stopRange = catima::range(proj, m_material); //get the total range for startEnergy -> 0, returns g/cm^2!
		proj.T = finalEnergy/proj.A;
		double finalRange = catima::range(proj, m_material); //get the range bound from the final energy;
		return (stopRange - finalRange) * densityInv * 0.01; //convert to m
	}

	//ZP, AP: projectile isotope, energy: MeV, pathLength: meters
	double Target::GetAngularStraggling(int zp, int ap, double energy, double pathLength)
	{
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = energy/proj.A;
		m_material.thickness_cm(pathLength * 100.0);
		proj.T = energy/proj.A;
		return catima::angular_straggling(proj, m_material);
	}

}