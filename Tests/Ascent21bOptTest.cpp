// vim:ts=2:et
//===========================================================================//
//                        "Tests/Ascent21bOptTest.cpp":                      //
//        Optimal Ascent of "Soyuz-2.1b" to Low Earth Orbit (via NOMAD)      //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent21b.h"
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>
#include <string>

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
// NB: As "AscentTest", plus optional "Section.Key=Value" args over-riding the
// Config.ini entries (eg "Mission.LaunchPad=Baykonur_31_6"):
//
int main(int argc, char* argv[])
{
  using namespace SpaceBallistics;
  using namespace std;

  if (argc < 2)
  {
    cerr << "PARAMETERS: Config.ini [Section.Key=Value ...]" << endl;
    return 1;
  }
  string configIni(argv[1]);

  try
  {
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(configIni, pt);

    for (int i = 2; i < argc; ++i)
    {
      string arg(argv[i]);
      size_t eq = arg.find('=');
      if (eq == string::npos || eq == 0)
      {
        cerr << "ERROR: Invalid Over-Ride: " << arg << endl;
        return 1;
      }
      pt.put(arg.substr(0, eq), arg.substr(eq + 1));
    }

    // NB: The log goes to "cout": not at its default 6 digits, which round
    // the Down-Range (~1800 km from the Insertion near the Lift-Off) to 10 m,
    // but as in "Ascent21bTest":
    cout.precision(15);

    auto [optRes, finalRunRes] =
      Ascent21b::FindOptimalAscentCtls(pt, &cout);

    if (!bool(optRes))
    {
      // Optimisation Failed:
      cout << "FIG!" << endl;
      return 1;
    }

    // If OK: Output the optimal params (or the least infeasible ones, if so
    // requested by "Technical.WithInfeasibleResult", or the fixed ones):
    string status(optRes.value().m_status);
    cout << ((status == "FEASIBLE")
             ? "OPTIMAL SOLUTION FOUND:"
             : (status == "NOT OPTIMISED")
             ? "FIXED PARAMS (NOT OPTIMISED):"
             : "NO FEASIBLE SOLUTION: LEAST INFEASIBLE ONE (NOT CLOSED):")
         << endl;
    cout << optRes.value()                  << endl;

    if (bool(finalRunRes))
    {
      auto const& frr = finalRunRes.value();
      // The closure (as in "Ascent21b::FindBurnTime3"): the Lift-Off State at
      // the Lift-Off time:
      Time dT     = frr.m_T + optRes.value().m_Tflight;
      bool closed =
        frr.m_rc != Ascent21b::RunRC::Error                  &&
        frr.m_VT <= VelK(pt.get<double>("Opt.MaxStartV"))     &&
        Abs(dT)  <= Time(pt.get<double>("Opt.MaxEndTimeMiss"));
      cout
        << (closed ? "FINAL RUN (CLOSED):" : "FINAL RUN (NOT CLOSED):")
        << endl
        << "RESULT    = " << Ascent21b::Base::ToString(frr.m_rc) << endl
        << "startTime = " << frr.m_T        << endl
        << "endTimeMiss = " << dT          << endl
        << "startL    = " << frr.m_LT       << endl
        << "startV    = " << frr.m_VT       << endl
        << "startMass = " << frr.m_mT       << endl
        << "maxQ      = " << frr.m_maxQ     << endl
        << "sepQ      = " << frr.m_sepQ     << endl
        << "maxLongG  = " << frr.m_maxLongG << endl
        << "payLoad   = " << optRes.value().m_payLoadMass        << endl
        << "T3        = " << optRes.value().m_T3                 << endl
        << "reserve3  = " << optRes.value().m_reserve3           << endl
        << endl;
    }
  }
  catch (exception const& exn)
  {
    cerr << "ERROR: " << exn.what() << endl;
    return 1;
  }
  return  0;
}
