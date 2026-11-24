// vim:ts=2:et
//===========================================================================//
//                         "Tests/Ascent21bTest.cpp":                        //
//   Bwd Ascent of "Soyuz-2.1b": Consistency with the Stage Models, Closure  //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent21b.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

namespace
{
  using namespace SpaceBallistics;
  namespace EAM = EarthAtmosphereModel;
  using     B   = Soyuz21b_Stage1_Booster<'B'>;

  constexpr FlightTime LiftOffFT(0.0_sec, TT(0.0_sec));

  //=========================================================================//
  // A Row of the Level-4 Log ("Ascent2" columns 1..17, plus (18) Spent):    //
  //=========================================================================//
  struct Row
  {
    double      m_t, m_h, m_mass, m_thrust, m_burnRate, m_p, m_spent;
    std::string m_mode;
  };

  //=========================================================================//
  // Direct Sum of the Stage Models (the Head Unit and PayLoad excluded):    //
  //=========================================================================//
  // A Stage is attached until (but not including) its Separation Time:
  //
  std::tuple<Mass, MassRate, Force> DirectSum(Time a_tau, Pressure a_p)
  {
    FlightTime ft = LiftOffFT + a_tau;
    Mass       m;
    MassRate   mDot;
    Force      f;
    auto add =
      [&m, &mDot, &f](StageDynParams<LVSC::Soyuz21b> const& a_dp)
      { m += a_dp.m_fullMass; mDot += a_dp.m_fullMassDot;
        f += a_dp.m_thrust[0]; };

    Density rho = EAM::Rho0;   // Not used with 0 velocity
    if (a_tau < Time(SC::Stage1SepTime))
    {
      constexpr B::VernDeflections vd0;
      constexpr ME::VelVE          v0;
      constexpr Angle_deg          fin0;
      add(Soyuz21b_BlockB::GetDynParams(ft, vd0, a_p, rho, v0, v0, fin0));
      add(Soyuz21b_BlockV::GetDynParams(ft, vd0, a_p, rho, v0, v0, fin0));
      add(Soyuz21b_BlockG::GetDynParams(ft, vd0, a_p, rho, v0, v0, fin0));
      add(Soyuz21b_BlockD::GetDynParams(ft, vd0, a_p, rho, v0, v0, fin0));
    }
    if (a_tau < Time(SC::Stage2SepTime))
    {
      constexpr S2::VernDeflections vd0;
      add(S2::GetDynParams(ft, vd0, a_p));
    }
    S3::ChamberDeflections const cd0 {};
    add(S3::GetDynParams(ft, cd0));
    return std::make_tuple(m, mDot, f);
  }

  int NFailed = 0;

  void Check(bool a_ok, std::string const& a_what)
  {
    std::cout << (a_ok ? "PASS: " : "FAIL: ") << a_what << std::endl;
    if (!a_ok)
      ++NFailed;
  }

  double RelErr(double a_x, double a_y)
    { return std::fabs(a_x - a_y) / std::max(std::fabs(a_y), 1.0); }

  // NB: From "Src/Missions/Ascent2.ini" (the "Ascent2" Ctl limits and the
  // Lift-Off Equivalent Velocity constraint):
  constexpr Angle_deg MaxAoA2   = Angle_deg(20.0);
  constexpr Angle_deg MaxAoA1   = Angle_deg(2.0);
  constexpr VelK      MaxStartV = VelK(0.03);

  constexpr Mass      PayLoad   = 8000.0_kg;
  // XXX: Not the 0.01 sec of "Ascent2.ini": the Stage Models assert (with
  // asserts on: all but the UnChecked mode) in "ToricSegm::PropMoICompsUp"
  // (rounding in the MoI of a nearly-empty toric tank segment) within ~1 msec
  // windows around Flight Times 42.646, 54.958, 68.778 (Boosters, Stage2) and
  // 150.321 sec (Stage2).
  // The RKF5 nodes are anchored at the Segment tops (Stage Model events), so
  // they are the same in every run;  with 0.025 sec, none falls into those
  // windows (with 0.01 sec, eg 54.9575 does):
  constexpr Time      Step      = 0.025_sec;

  // The Dim-Less AoA params (aAoAHat2, bAoAHat2, aAoAHat1, bAoAHat1): the
  // "Ascent2" defaults (all 0, ie a zero-AoA gravity turn), and a set which
  // brings the Plesetsk case close to the Lift-Off state. XXX: The latter is
  // hand-picked by a coarse grid scan on that case; it is NOT a fit (which is
  // to be done by NOMAD, as in "Ascent2-Opt"):
  struct CtlSet
  {
    char const* m_name;
    double      m_hats[4];
  };
  constexpr CtlSet CtlSets[]
  {
    { "Ascent2 defaults", { 0.0,   0.0,   0.0, 0.0 } },
    { "Hand-picked",      { 0.425, 0.295, 0.0, 0.0 } }
  };

  void PrintRes(Ascent21b::CloseRes const& a_cr)
  {
    using     std::cout;
    using     std::endl;
    auto const& r = a_cr.m_res;
    cout << "#   Verdict   = " << (a_cr.m_closed ? "CLOSED" : "NOT CLOSED")
         << endl
         << "#   RC        = " << Ascent21b::Base::ToString(r.m_rc)    << endl
         << "#   T3        = " << a_cr.m_T3                      << endl
         << "#   Tflight   = " << - r.m_T + a_cr.m_dT            << endl
         << "#   dT        = " << a_cr.m_dT
         << "  (end time - lift-off time)"                       << endl
         << "#   VT        = " << r.m_VT
         << "  (Equiv. Velocity @ H=0; limit " << MaxStartV << ")" << endl
         << "#   LT        = " << r.m_LT                         << endl
         << "#   mT        = " << r.m_mT                         << endl
         << "#   reserve3  = " << a_cr.m_reserve3                << endl
         << "#   maxQ      = " << r.m_maxQ                       << endl
         << "#   sepQ      = " << r.m_sepQ                       << endl
         << "#   maxLongG  = " << r.m_maxLongG                   << endl
         << "#   nRuns     = " << a_cr.m_nRuns                   << endl;
  }
}

//===========================================================================//
// "main":                                                                   //
//===========================================================================//
int main()
{
  using namespace std;

  //-------------------------------------------------------------------------//
  // The Plesetsk Case (Pad 43/4): Closure with the full Log:                //
  //-------------------------------------------------------------------------//
  ostringstream log;
  log.precision(15);

  Ascent21b proto
    (SC::Stage3BurnDur, MaxAoA2, MaxAoA1, PayLoad, 290.0_km, 340.0_km,
     Angle_deg(64.0), Pad_Plesetsk_43_4, Step, &log, 4);
  {
    double const* h = CtlSets[1].m_hats;
    proto.SetCtlParams(h[0], h[1], h[2], h[3]);
  }

  Ascent21b::CloseRes cr = proto.FindBurnTime3(MaxStartV);
  Check(cr.m_res.m_rc != Ascent21b::RunRC::Error,
        string("Plesetsk_43_4: closure search ran, RC=") +
        Ascent21b::Base::ToString(cr.m_res.m_rc));

  Time Tflight =
    Time(SC::Stage2CutOffTime) - Time(SC::LiftOffTime) + cr.m_T3;

  //-------------------------------------------------------------------------//
  // Parse the Log (of the final run: the search runs are at level 2):       //
  //-------------------------------------------------------------------------//
  vector<Row>                    rows;
  vector<pair<double, string>>   events;
  {
    istringstream is(log.str());
    string        line;
    while (getline(is, line))
    {
      if (line.rfind("# t=", 0) == 0)
      {
        // Event: "# t=<t> sec: <Name>: ..."
        size_t sec = line.find(" sec: ");
        size_t end = line.find(": ", sec + 6);
        if (sec != string::npos && end != string::npos)
          events.emplace_back(stod(line.substr(4, sec - 4)),
                              line.substr(sec + 6, end - sec - 6));
        continue;
      }
      if (line.empty() || line[0] == '#')
        continue;
      istringstream ls(line);
      Row    r;
      double L, Vr, Vhor, V, psi, aoa, Q, M, longG, dt;
      ls >> r.m_t >> r.m_h >> L >> Vr >> Vhor >> V >> psi >> aoa >> r.m_mass
         >> r.m_mode >> r.m_thrust >> r.m_burnRate >> Q >> M >> longG
         >> r.m_p >> dt >> r.m_spent;
      if (ls)
        rows.push_back(r);
    }
  }
  Check(rows.size() > 1000, "Level-4 log parsed: " + to_string(rows.size()) +
        " rows of 18 columns");
  if (rows.size() <= 1000)
  {
    cout << log.str();
    return 1;
  }

  auto eventTime =
    [&events](string const& a_name) -> double
    {
      for (auto const& ev: events)
        if (ev.second == a_name)
          return ev.first;
      return NaN<double>;
    };

  //-------------------------------------------------------------------------//
  // Staging at the Stage Models' Times:                                     //
  //-------------------------------------------------------------------------//
  struct { char const* m_name; Time m_tau; } const expEvents[]
  {
    { "Stage1 Full Thrust",              Time(B ::FullThrustTime)   },
    { "Stage1 Throttle-Down",            Time(B ::ThrottlTime)      },
    { "Stage1 Cut-Off",                  Time(SC::Stage1CutOffTime) },
    { "Stage1 Separation",               Time(SC::Stage1SepTime)    },
    { "Stage2 Throttle-Down",            Time(S2::ThrottlTime)      },
    { "Stage2 Main Engine Cut-Off",      Time(S2::MainCutOffTime)   },
    { "Stage2 Cut-Off, Stage3 Ignition", Time(SC::Stage2CutOffTime) },
    { "Stage2 Separation",               Time(SC::Stage2SepTime)    },
    { "Stage3 Aft Section Jettison",     Time(SC::Stage3AftJetTime) },
    { "Stage3 Cut-Off, Orbital Insertion", Tflight                  }
  };
  vector<double> evTaus;
  bool reachedLiftOff = !std::isnan(eventTime("Lift-Off"));
  for (auto const& ee: expEvents)
  {
    double te  = eventTime(ee.m_name) + Tflight.Magnitude();
    double tm  = ee.m_tau.Magnitude();
    evTaus.push_back(tm);
    // Events below the end of the Bwd run are not reached, which is OK only
    // for the earliest ones:
    if (std::isnan(te))
    {
      cout << "# Not reached: " << ee.m_name << endl;
      continue;
    }
    Check(std::fabs(te - tm) < 1e-9, string(ee.m_name) + " @ tau=" +
          to_string(te) + " sec (Stage Models: " + to_string(tm) + " sec)");
  }
  double tauF = eventTime("Fairing Separation") + Tflight.Magnitude();
  Check(Time(SC::FairingJetTimeMin).Magnitude() <= tauF &&
        tauF <= Time(SC::FairingJetTimeMax).Magnitude(),
        "Fairing Separation @ tau=" + to_string(tauF) +
        " sec, within [FairingJetTimeMin, FairingJetTimeMax]");

  //-------------------------------------------------------------------------//
  // Every Row vs the Direct Sum of the Stage Models:                        //
  //-------------------------------------------------------------------------//
  // Mass, BurnRate and Thrust at the Row's Fwd Time and ambient pressure; the
  // Mode, and the Spent Propellant (the 4th State component) vs the Stage
  // Models' mass drop minus the dry masses jettisoned after the Row's time.
  // Rows at (or within rounding of) an event time are skipped:
  //
  constexpr S2::VernDeflections vd2;
  B::VernDeflections const      vd1 {};
  ME::VelVE const               v0  {};
  Mass jettB = 4.0 * B::GetDynParams
    (LiftOffFT + Time(SC::Stage1SepTime), vd1, EAM::P0, EAM::Rho0, v0, v0,
     Angle_deg(0.0)).m_fullMass;
  Mass jettC = S2::GetDynParams
    (LiftOffFT + Time(SC::Stage2SepTime), vd2, Pressure(0.0)).m_fullMass;
  Mass jettA = S3::AftMass;
  Mass jettF = SH::FairingLargeMass;
  Mass head  = SH::InterStageLargeMass + PayLoad;

  double mIns   = rows.front().m_mass;
  double maxEm  = 0.0, maxEf = 0.0, maxEmu = 0.0, maxEs = 0.0, tauEs = 0.0;
  bool   modeOK = true;
  int    nCmp   = 0;

  for (Row const& r: rows)
  {
    double tau = r.m_t + Tflight.Magnitude();
    bool   skip = (r.m_t == 0.0 || tau < 1e-9);
    for (double te: evTaus)
      skip = skip || std::fabs(tau - te) < 1e-6;
    skip   = skip || std::fabs(tau - tauF) < 1e-6;
    if (skip)
      continue;
    ++nCmp;

    Pressure p   = std::min(Pressure(r.m_p), EAM::P0);
    auto [ms, mDot, fs] = DirectSum(Time(tau), p);
    bool     fair = r.m_t < tauF - Tflight.Magnitude();
    Mass     m    = ms + head + (fair ? jettF : 0.0_kg);
    double   fkg  = (fs / g0).Magnitude();

    maxEm  = std::max(maxEm,  RelErr(r.m_mass,     m.Magnitude()));
    maxEf  = std::max(maxEf,  RelErr(r.m_thrust,   fkg));
    maxEmu = std::max(maxEmu, RelErr(r.m_burnRate, - mDot.Magnitude()));

    // Expected Mode:
    char const* em =
      (tau > Time(SC::Stage2SepTime)   .Magnitude()) ? "Burn3"  :
      (tau > Time(SC::Stage2CutOffTime).Magnitude()) ? "Burn23" :
      (tau > Time(SC::Stage1SepTime)   .Magnitude()) ? "Burn2"  : "Burn12";
    if (r.m_mode != em)
      modeOK = false;

    // Spent Propellant Book-Keeping:
    Mass jett =
      (tau < Time(SC::Stage3AftJetTime).Magnitude() ? jettA : 0.0_kg) +
      (tau < Time(SC::Stage2SepTime)   .Magnitude() ? jettC : 0.0_kg) +
      (tau < Time(SC::Stage1SepTime)   .Magnitude() ? jettB : 0.0_kg) +
      (fair                                         ? jettF : 0.0_kg);
    double sp = r.m_mass - mIns - jett.Magnitude();
    if (std::fabs(sp - r.m_spent) > maxEs)
    {
      maxEs = std::fabs(sp - r.m_spent);
      tauEs = tau;
    }
  }
  {
    Mass m0 = get<0>(DirectSum(Tflight, Pressure(0.0))) + head;
    cout << "# Insertion Mass: Log=" << mIns << " kg, StageModels+Head+"
            "PayLoad=" << m0.Magnitude() << " kg" << endl;
    Check(RelErr(mIns, m0.Magnitude()) < 1e-9,
          "Insertion mass = sum of Stage Models + InterStage + PayLoad");
  }
  Check(nCmp > 1000 && maxEm < 1e-12 && maxEf < 1e-9 && maxEmu < 1e-9,
        to_string(nCmp) + " rows: mass, burn rate, thrust = sum of Stage "
        "Models (max rel errs " + to_string(maxEm) + ", " + to_string(maxEmu)
        + ", " + to_string(maxEf) + ")");
  Check(modeOK, "Flight Modes follow the Stage Models' separation times");
  Check(maxEs < 1e-3,
        "Spent propellant = Stage Models' mass drop - jettisoned (max err " +
        to_string(maxEs) + " kg @ tau=" + to_string(tauEs) + " sec)");

  //-------------------------------------------------------------------------//
  // The Closure:                                                            //
  //-------------------------------------------------------------------------//
  cout << "# Plesetsk_43_4, i=64 deg, 290 x 340 km, PayLoad=" << PayLoad
       << ", Ctls: " << CtlSets[1].m_name << ": reached Lift-Off time: "
       << (reachedLiftOff ? "yes" : "no") << endl;
  PrintRes(cr);
  if (cr.m_res.m_rc != Ascent21b::RunRC::Error)
  {
    Mass m0 = get<0>(DirectSum(0.0_sec, EAM::P0)) + head + jettF;
    cout << "#   (Lift-Off mass of the Stage Models + Head + PayLoad = "
         << m0 << ")" << endl;
  }

  //-------------------------------------------------------------------------//
  // All the library's "Soyuz-2" Pads (LaunchPads.hpp):                      //
  //-------------------------------------------------------------------------//
  struct PadCase
  {
    char const*             m_name;
    Soyuz2_LaunchPad const* m_pad;
    double                  m_incl;   // deg
  };
  PadCase const pads[]
  {
    { "Plesetsk_43_3",  &Pad_Plesetsk_43_3,  64.0 },
    { "Plesetsk_43_4",  &Pad_Plesetsk_43_4,  64.0 },
    { "Baykonur_31_6",  &Pad_Baykonur_31_6,  51.6 },
    { "Vostochny_1S",   &Pad_Vostochny_1S,   97.4 }
  };
  for (CtlSet  const& cs: CtlSets)
  for (PadCase const& pc: pads)
  {
    Ascent21b a
      (SC::Stage3BurnDur, MaxAoA2, MaxAoA1, PayLoad, 290.0_km, 340.0_km,
       Angle_deg(pc.m_incl), *pc.m_pad, Step, nullptr, 0);
    a.SetCtlParams(cs.m_hats[0], cs.m_hats[1], cs.m_hats[2], cs.m_hats[3]);
    Ascent21b::CloseRes c = a.FindBurnTime3(MaxStartV);
    cout << "# " << pc.m_name << " (" << cs.m_name << "): i=" << pc.m_incl
         << " deg, Elev=" << pc.m_pad->Elevation().Magnitude() << " m: "
         << (c.m_closed ? "CLOSED" : "NOT CLOSED")
         << ", RC="       << Ascent21b::Base::ToString(c.m_res.m_rc)
         << ", T3="       << c.m_T3.Magnitude()
         << " sec, dT="   << c.m_dT.Magnitude()
         << " sec, VT="   << c.m_res.m_VT.Magnitude()
         << " km/sec, reserve3=" << c.m_reserve3.Magnitude()
         << " kg, maxQ="  << c.m_res.m_maxQ.Magnitude()
         << " Pa, runs="  << c.m_nRuns << endl;
    Check(c.m_res.m_rc != Ascent21b::RunRC::Error,
          string(pc.m_name) + ": closure search ran without Error");
  }

  cout << (NFailed == 0 ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED") << endl;
  return NFailed == 0 ? 0 : 1;
}
