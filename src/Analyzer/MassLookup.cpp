#include "MassLookup.h"

#include <fstream>
#include <sstream>

namespace AnasenAnalyzer {

    MassLookup* MassLookup::s_instance = new MassLookup();

	MassLookup::MassLookup()
	{
		std::ifstream massfile("etc/mass.txt");
		if(massfile.is_open())
		{
			std::string junk, element;
			double atomicMassBig, atomicMassSmall, isotopicMass;
			getline(massfile,junk);
			getline(massfile,junk);
			NucData data;
			while(massfile>>junk)
			{
				massfile>>data.Z>>data.A>>element>>atomicMassBig>>atomicMassSmall;
				data.mass = (atomicMassBig + atomicMassSmall*1e-6 - data.Z*s_electronMass)*s_u2MeV;
				data.isotopicSymbol = std::to_string(data.A) + element;
				m_dataMap[GetNucleusUUID(data.Z, data.A)] = data;
			}
		}
	}
	
	MassLookup::~MassLookup() {}
	
	//Returns nuclear mass in MeV
	double MassLookup::FindMass(uint32_t Z, uint32_t A)
	{
		auto data = m_dataMap.find(GetNucleusUUID(Z, A));
		if(data == m_dataMap.end())
			return s_invalidMass;
	
		return data->second.mass;
	}
	
	//returns element symbol
	std::string MassLookup::FindSymbol(uint32_t Z, uint32_t A)
	{
		auto data = m_dataMap.find(GetNucleusUUID(Z, A));
		if(data == m_dataMap.end())
			return s_invalidSymbol;
	
		return data->second.isotopicSymbol;
	}

	const NucData& MassLookup::FindData(uint32_t Z, uint32_t A)
	{
		auto data = m_dataMap.find(GetNucleusUUID(Z, A));
		if(data == m_dataMap.end())
			return m_invalidResult;
	
		return data->second;
	}
}