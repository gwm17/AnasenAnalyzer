#ifdef __CINT__
 
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;

#pragma link C++ struct AnasenAnalyzer::Silicon_Event+;
#pragma link C++ struct AnasenAnalyzer::PropCounter_Event+;
#pragma link C++ struct AnasenAnalyzer::CsI_Event+;
#pragma link C++ struct SiDetector+;
#pragma link C++ struct SiEvent+;
#pragma link C++ struct PCEvent+;
#pragma link C++ struct CsIEvent+;
#pragma link C++ struct CsIHit::SortByCsI+;
#pragma link C++ struct SiHit::SortByDetector+;
#pragma link C++ struct SiHit::SortByHit+;
#pragma link C++ struct PCHit::SortByPC+;

#pragma link C++ class std::vector<AnasenAnalyzer::Silicon_Event>+;
#pragma link C++ class std::vector<AnasenAnalyzer::PropCounter_Event>+;
#pragma link C++ class std::vector<AnasenAnalyzer::CsI_Event>+;
#pragma link C++ class std::vector<SiHit::SortByHit>+;
#pragma link C++ class std::vector<SiHit::SortByDetector>+;
#pragma link C++ class std::vector<PCHit::SortByPC>+;
#pragma link C++ class std::vector<CsIHit::SortByCsI>+;

#pragma link C++ defined_in "PCHit.h";
#pragma link C++ defined_in "SiHit.h";
#pragma link C++ defined_in "CsIHit.h";

#endif