#ifndef MASS_LOOKUP_H
#define MASS_LOOKUP_H

#include <unordered_map>
#include <cstdint>
#include <string>

namespace AnasenAnalyzer {

    static constexpr uint32_t GetNucleusUUID(uint32_t z, uint32_t a)
    {
        return z >= a ? z*z + z + a : a*a + z;
    }
	struct NucData
    {
        double mass = -1.0; //MeV
        uint32_t Z = 0;
        uint32_t A = 0;
        std::string isotopicSymbol = "Invalid";
    };

	class MassLookup
	{
	public:
		~MassLookup();
		double FindMass(uint32_t Z, uint32_t A);
		double FindMassU(uint32_t Z, uint32_t A) { return FindMass(Z, A) / s_u2MeV; }
		std::string FindSymbol(uint32_t Z, uint32_t A);
		const NucData& FindData(uint32_t Z, uint32_t A);

	    static bool IsInvalidSymbol(const std::string& symbol) { return symbol == s_invalidSymbol; }
		static bool IsInvalidMass(double mass) { return mass == s_invalidMass; }
		static MassLookup& GetInstance() { return *s_instance; }
		
	
	private:
		MassLookup();

		static MassLookup* s_instance;
		std::unordered_map<uint32_t, NucData> m_dataMap;
		NucData m_invalidResult;
	
		//constants
		static constexpr double s_u2MeV = 931.4940954;
		static constexpr double s_electronMass = 0.000548579909; //MeV
        static constexpr double s_invalidMass = -1.0;
        static constexpr char s_invalidSymbol[] = "Invalid";
	    
	};
}

#endif