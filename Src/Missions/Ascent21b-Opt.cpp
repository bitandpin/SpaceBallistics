// vim:ts=2:et
//===========================================================================//
//                      "Src/Missions/Ascent21b-Opt.cpp":                    //
//  Ascent-to-Orbit of the "Soyuz-2.1b" LV: Parametric Optimisation (NOMAD)  //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent21b.h"
#include "SpaceBallistics/Missions/MkNOMADParams.hpp"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.h"
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <cstdlib>

namespace SpaceBallistics
{
//===========================================================================//
// "ModifyLVParams":                                                         //
//===========================================================================//
// This method allows us to modify consistently the Mass and Thrust params in
// the "Ascent21b" object. Intended for use with Optimisation Procedures:
//
void Ascent21b::ModifyLVParams
(
  // Mass Params:
  Mass      a_payLoadMass,

  // Ctl Params:
  double    a_aAoAHat2,   double a_bAoAHat2,
  Time      a_T3,
  double    a_aAoAHat1,   double a_bAoAHat1
)
{
  // Verify the Mass Params (the Ctl Params are verified later in "SetBurn-
  // Time3" and "SetCtlParams"):
  if (IsNeg(a_payLoadMass))
    throw std::invalid_argument
          ("Ascent21b::ModifyLVParams: Invalid Mass Param");

  //-------------------------------------------------------------------------//
  // Set the Mass Params:                                                    //
  //-------------------------------------------------------------------------//
  // Over-All:
  m_payLoadMass  = a_payLoadMass;

  // Soyuz-2.1b: the cached Stage Models' sum includes the PayLoad, so inval-
  // idate it:
  m_stack.m_tau  = Time(NAN);

  // The Ctls:
  SetBurnTime3(a_T3);
  SetCtlParams(a_aAoAHat2, a_bAoAHat2, a_aAoAHat1, a_bAoAHat1);
}

//===========================================================================//
// "OptRes" Non-Default Ctor:                                                //
//===========================================================================//
Ascent21b::OptRes::OptRes(Ascent21b const& a_asc, char const* a_status)
:
  // For Stage2:
  m_aAoAHat2   (a_asc.m_aAoAHat2),
  m_bAoAHat2   (a_asc.m_bAoAHat2),
  m_aAoA2      (a_asc.m_aAoA2),
  m_bAoA2      (a_asc.m_bAoA2),

  // For Stage1:
  m_aAoAHat1   (a_asc.m_aAoAHat1),
  m_bAoAHat1   (a_asc.m_bAoAHat1),
  m_aAoA1      (a_asc.m_aAoA1),
  m_bAoA1      (a_asc.m_bAoA1),

  // Over-All:
  m_payLoadMass(a_asc.m_payLoadMass),

  // Soyuz-2.1b: the AoA Domains and the Stage3 BurnTime:
  m_T2         (a_asc.m_T2),
  m_T1         (a_asc.m_T1),
  m_T3         (a_asc.m_T3),
  m_Tflight    (a_asc.m_Tflight),
  m_reserve3   (a_asc.GetReserve3()),
  m_status     (a_status)
{}

//===========================================================================//
// "OptRes" Output:                                                          //
//===========================================================================//
std::ostream& operator<< (std::ostream& a_os, Ascent21b::OptRes const& a_res)
{
  return
  a_os
    << "\tstatus         = "     << a_res.m_status
    // Dim-Less Params (NP):
    << "\n\taAoAHat2       = "   << a_res.m_aAoAHat2
    << "\n\tbAoAHat2       = "   << a_res.m_bAoAHat2
    << "\n\taAoAHat1       = "   << a_res.m_aAoAHat1
    << "\n\tbAoAHat1       = "   << a_res.m_bAoAHat1
    // Dimensioned Params:
    << "\n\tT2             = "   << a_res.m_T2
    << "\n\taAoA2          = "   << a_res.m_aAoA2
    << "\n\tbAoA2          = "   << a_res.m_bAoA2
    << "\n\tT3             = "   << a_res.m_T3
    << "\n\tTflight        = "   << a_res.m_Tflight
    << "\n\treserve3       = "   << a_res.m_reserve3
    << "\n\tT1             = "   << a_res.m_T1
    << "\n\taAoA1          = "   << a_res.m_aAoA1
    << "\n\tbAoA1          = "   << a_res.m_bAoA1
    << "\n\tpayLoadMass    = "   << a_res.m_payLoadMass
    // XXX: Ctl Eqs (the AoA laws only):
    << "\n\t# t2  = -"           << a_res.m_T2.Magnitude()
    << " .. 0"
       "\n\tAoA2 := t2 * ("      << a_res.m_aAoA2.Magnitude()
    << " * t2 - ("               << a_res.m_bAoA2.Magnitude()      << "));"
    << "\n\t# tau1 = 0 .. "      << a_res.m_T1.Magnitude()
    << "\n\tAoA1 := tau1 * ("    << a_res.m_aAoA1.Magnitude()
    << " * tau1 + ("             << a_res.m_bAoA1.Magnitude()      << ")); "
    << std::endl;
}

//===========================================================================//
// "Ascent21b::NOMADEvaluator": Helper Class used in NOMAD Optimisation:     //
//===========================================================================//
class Ascent21b::NOMADEvaluator final: public NOMAD::Evaluator
{
private:
  //=========================================================================//
  // Data Flds:                                                              //
  //=========================================================================//
  Ascent21b const*    const  m_proto;
  std::array<bool,NP> const  m_actOpts;    // Flags for Active Opt Params
  Mass                const  m_massScale;  // For "payLoadMassRel"
  // Optimisation Limits (Constraints):
  VelK                const  m_startVLimit;
  Pressure            const  m_QLimit;     // NAN -> disabled
  Pressure            const  m_sepQLimit;  // ditto
  double              const  m_longGLimit; // ditto
  Time               const  m_endTimeTol; // On |End Time Miss|

public:
  //=========================================================================//
  // Non-Default Ctor, Dtor:                                                 //
  //=========================================================================//
  NOMADEvaluator
  (
    std::shared_ptr<NOMAD::EvalParameters> const& a_params,
    Ascent21b           const*                    a_proto,
    // Flags indicating which variables are optimisation args:
    std::array<bool,NP> const&                    a_act_opts,
    Mass                                          a_mass_scale,
    // Optimisation Limits (Constraints):
    VelK                                          a_startV_limit,
    Pressure                                      a_Q_limit,
    Pressure                                      a_sepQ_limit,
    double                                        a_longG_limit,
    Time                                          a_end_time_tol
  )
  : NOMAD::Evaluator(a_params, NOMAD::EvalType::BB),
    m_proto         (a_proto),
    m_actOpts       (a_act_opts),
    m_massScale     (a_mass_scale),
    m_startVLimit   (a_startV_limit),
    m_QLimit        (a_Q_limit),
    m_sepQLimit     (a_sepQ_limit),
    m_longGLimit    (a_longG_limit),
    m_endTimeTol    (a_end_time_tol)
  {
    assert(m_proto != nullptr);

    if (!(IsPos(m_massScale)   && IsPos(m_startVLimit) && IsPos(m_QLimit) &&
          IsPos(m_sepQLimit)   && IsPos(m_longGLimit)  &&
          IsPos(m_endTimeTol)))
      throw std::invalid_argument
            ("Ascent21b::NOMADEvaluator::Ctor: Invalid Limit(s)");
  }

  ~NOMADEvaluator() override {}

  //=========================================================================//
  // "eval_x": The Actual Evaluation Method for NOMAD:                       //
  //=========================================================================//
  bool eval_x
  (
    NOMAD::EvalPoint&    a_x,  // Not "const" -- the result is also set here
    NOMAD::Double const& UNUSED_PARAM(a_hMax),
    bool&                a_countEval
  )
  const override
  {
    //-----------------------------------------------------------------------//
    // Defaults for all Optimisation Params (unless over-written):           //
    //-----------------------------------------------------------------------//
    // The over-all list of params is:
    // [aAoAHat2, bAoAHat2, T3Rel, aAoAHat1, bAoAHat1, payLoadMassRel],
    // where all params are Dim-Less:
    //
    double aAoAHat2     = m_proto->m_aAoAHat2;
    double bAoAHat2     = m_proto->m_bAoAHat2;
    Time   T3           = m_proto->m_T3;
    double aAoAHat1     = m_proto->m_aAoAHat1;
    double bAoAHat1     = m_proto->m_bAoAHat1;
    Mass   payLoadMass  = m_proto->m_payLoadMass;

    //-----------------------------------------------------------------------//
    // Extract the actual params from "a_x":                                 //
    //-----------------------------------------------------------------------//
    int j = 0;
    for (int i = 0; i < NP; ++i)
    {
      if (!m_actOpts[size_t(i)])
        continue;

      // Otherwise: the "i"th over-all idx is mapped to the "j"th idx in "a_x":
      double xj = a_x[size_t(j)].todouble();

      switch (i)
      {
      case  0: aAoAHat2    = xj; break;
      case  1: bAoAHat2    = xj; break;
      case  2: T3          = xj * SC::Stage3BurnDur; break;
      case  3: aAoAHat1    = xj; break;
      case  4: bAoAHat1    = xj; break;
      case  5: payLoadMass = xj * m_massScale;       break;
      default: assert(false);
      }
      ++j;
    }
    //-----------------------------------------------------------------------//
    // For thread-safety, construct a new "Ascent21b" obj from "m_proto":    //
    //-----------------------------------------------------------------------//
    // NB: This is a MUST: "Ascent21b" has a "mutable" Stage Models' cache:
    Ascent21b asc(*m_proto);
    RunRes    res;

    // Soyuz-2.1b: an invalid "T3" or an exceeded MaxAoA throws an exception
    // which is not handled inside "Run", so it is caught here, as a failed
    // evaluation:
    try
    {
      // Install the possibly-modified params in "asc":
      asc.ModifyLVParams
      (
        // Mass Params:
        payLoadMass,
        // Ctl Params:
        aAoAHat2,    bAoAHat2,  T3,
        aAoAHat1,    bAoAHat1
      );

      //---------------------------------------------------------------------//
      // Run the Integrator!                                                 //
      //---------------------------------------------------------------------//
      res = asc.Run();
    }
    catch (std::exception const& exn)
    {
      if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 1)
#       pragma omp critical(Output)
        *(m_proto->m_os) << "# Ascent21b::NOMADEvaluator: ERROR: "
                         << exn.what() << std::endl;
      return false;
    }

    // IMPORTANT: NOMAD allows us to indicate that evaluation has failed:
    if (res.m_rc == RunRC::Error)
      return false;

    //-----------------------------------------------------------------------//
    // Push the results back to NOMAD in string form:                        //
    //-----------------------------------------------------------------------//
    char  buff[512];
    char* curr = buff;

    // XXX: The Objective Function Value: StartMass or (-PayLoadMass), both to
    // be minimised. Which one is to be used, depends on whether the PayLoad-
    // Mass is the list of Optimisation Params. Soyuz-2.1b: the StartMass is
    // fixed by the Stage Models, so the Stage3 reserve is subtracted from it
    // (ie the StartMass with the unused Stage3 Propellant off-loaded):
    bool objIsStartMass = !m_actOpts[NP-1];
    Mass objMass        =
      objIsStartMass
      ? (res.m_mT - asc.GetReserve3())
      : (- asc.m_payLoadMass);
    curr += sprintf(curr, "%.16e", objMass.Magnitude());

    // Constraints:
    // StartV (or a StartV equivalent taking StartH into account):
    VelK VT = std::max(res.m_VT, VelK(0.0));
    assert(!(IsNeg(res.m_maxQ) || IsNeg(res.m_sepQ) || IsNeg(res.m_maxLongG)));

    // Soyuz-2.1b: the End Time Miss: the Bwd integration must end at the
    // Lift-Off time (otherwise the Stage Models' masses at its end are not
    // the Lift-Off ones), as in "FindBurnTime3":
    Time dT = res.m_T + asc.m_Tflight;
    if (!(IsFinite(VT) && IsFinite(dT)))
      return false;

    curr += sprintf(curr, " %.16e", (VT - m_startVLimit).Magnitude());
    curr += sprintf(curr, " %.16e", (res.m_maxQ     - m_QLimit).Magnitude());
    curr += sprintf(curr, " %.16e", (res.m_sepQ - m_sepQLimit).Magnitude());
    curr += sprintf(curr, " %.16e", (res.m_maxLongG - m_longGLimit));
    curr += sprintf(curr, " %.16e", (Abs(dT) - m_endTimeTol).Magnitude());
    // Done!
    assert(size_t(curr - buff) <= sizeof(buff));

    if (m_proto->m_os != nullptr && m_proto->m_logLevel >= 2)
    {
#     pragma omp critical(Output)
      *(m_proto->m_os) << buff << std::endl;
    }
    // Set the results back in "a_x":
    a_x.setBBO(buff);
    a_countEval = true;
    return        true;
  }
};

//===========================================================================//
// "FindOptimalAscentCtls"                                                   //
//===========================================================================//
std::pair<std::optional<Ascent21b::OptRes>,
          std::optional<Ascent21b::Base::RunRes>> // If FinalRun is done
Ascent21b::FindOptimalAscentCtls
(
  boost::property_tree::ptree const& a_config,
  std::ostream*                      a_os       // May be NULL
)
{
  //-------------------------------------------------------------------------//
  // Bounds and Initial Vals for all Params:                                 //
  //-------------------------------------------------------------------------//
  // [aAoAHat2, bAoAHat2, T3Rel, aAoAHat1, bAoAHat1, payLoadMassRel]:
  // Soyuz-2.1b: of the "Ascent2" params, only the AoA Ctls and the PayLoad-
  // Mass are left: the Stage Models fix the thrusts, burn rates and masses,
  // and there is no Ballistic Gap. Instead, the Stage3 BurnTime (T3, the
  // closing param) is an Optimisation Param, and the closure (the Bwd integ-
  // ration ending at the Lift-Off time, within "Opt.MaxEndTimeMiss") is an
  // extra Constraint;  it replaces the search of "FindBurnTime3" which would
  // cost ~20 runs per evaluation. T3Rel is relative to the full burn of the
  // Stage Model, and payLoadMassRel to the Stage Models' Lift-Off Mass (incl
  // the Head Unit, but w/o the PayLoad, see below):
  double LoBounds[NP]
    { 0.0,      0.0,      0.5,   0.0,      0.0,      0.02 };

  double UpBounds[NP]
    { 1.0,      1.0,      1.0,   1.0,      1.0,      0.06 };

  double InitVals[NP]
    { 0.5,      0.5,      0.98,  0.5,      0.5,      0.04 };

  //-------------------------------------------------------------------------//
  // Open and Parse the Config.ini File:                                     //
  //-------------------------------------------------------------------------//
  // NB: Already parsed (see the overload below):
  boost::property_tree::ptree const& pt = a_config;

  // Get the LV and Mission params for the "prototype" "Ascent21b" obj:
  // Soyuz-2.1b: only the AoA limits, the rest come from the Stage Models:
  Angle_deg maxAoA2           (pt.get<double>("LV.MaxAoA2"));
  Angle_deg maxAoA1           (pt.get<double>("LV.MaxAoA1"));

  // Mission:
  LenK      perigee           (pt.get<double>("Mission.Perigee"));
  LenK      apogee            (pt.get<double>("Mission.Apogee" ));
  Angle_deg incl              (pt.get<double>("Mission.Inclination"));
  // Soyuz-2.1b: the Launch Site is a "Soyuz2_LaunchPad", not a Latitude:
  std::string padName       =  pt.get<std::string>("Mission.LaunchPad");

  Soyuz2_LaunchPad const* pad =
    (padName == "Plesetsk_43_3") ? &Pad_Plesetsk_43_3 :
    (padName == "Plesetsk_43_4") ? &Pad_Plesetsk_43_4 :
    (padName == "Baykonur_31_6") ? &Pad_Baykonur_31_6 :
    (padName == "Vostochny_1S" ) ? &Pad_Vostochny_1S  : nullptr;
  if (pad == nullptr)
    throw std::invalid_argument
          ("Ascent21b::FindOptimalAscentCtls: Invalid LaunchPad: " + padName);

  //-------------------------------------------------------------------------//
  // Soyuz-2.1b: Technical Params and the "prototype" obj First:             //
  //-------------------------------------------------------------------------//
  // The PayLoadMass scale (see below) comes from the "prototype" obj, so the
  // latter is created before the "Opt" params are read:
  //
  // Generic Params:
  Time   odeIntegrStep     (pt.get<double>("Technical.ODEIntegrStep"));
  bool   withFinalRun     = pt.get<bool>  ("Technical.WithFinalRun");
  int    finalRunLogLevel = pt.get<int>   ("Technical.FinalRunLogLevel");
  int    optMaxEvals      = pt.get<int>   ("Technical.OptMaxEvals");
  int    optLogLevel      = pt.get<int>   ("Technical.OptLogLevel");

  // NOMAD-Specific Params:
  int    optSeed          = pt.get<int>   ("Technical.NOMADSeed");
  bool   stopIfFeasible   = pt.get<bool>  ("Technical.NOMADStopIfFeasible");
  double useVNS           = pt.get<double>("Technical.NOMADUseVNS");
  if (useVNS < 0.0 || useVNS >= 1.0)      // 0: VNS not used
    throw std::invalid_argument ("NOMADUseVNS: The arg must be in [0..1)");
  bool   useMT            = pt.get<bool>  ("Technical.NOMADUseMT");

  // If no feasible point is found, use the least infeasible one (default: no
  // result):
  bool   withInfeasible   =
    pt.get<bool>("Technical.WithInfeasibleResult", false);

  // Create the "prototype" "Ascent21b" obj, with a zero PayLoad and the full
  // Stage3 burn for the moment; the actual vals are installed by "Modify-
  // LVParams" below:
  Ascent21b proto
  (
    SC::Stage3BurnDur,  maxAoA2, maxAoA1, 0.0_kg,
    perigee, apogee,    incl,    *pad,
    odeIntegrStep,      a_os,    optLogLevel
  );

  // The Mass Scale for "payLoadMassRel": the Stage Models' Lift-Off Mass,
  // incl the InterStage and the Fairing (the PayLoad is 0 in "proto"):
  Mass massScale =
    proto.GetStack(0.0_sec, EAM::P0, EAM::Rho0).m_mass + SH::FairingLargeMass;

  //-------------------------------------------------------------------------//
  // Now  get the "Opt" params:                                              //
  //-------------------------------------------------------------------------//
  // (*) If a param is give a numerical value in Config.ini, it is considered
  //     to be "fixed" (ie NOT an optimisation argument), and the corresp bit
  //     in "actOpts" is set to "false";
  // (*) if a param has the special value "OPTIMISE", is will serve as an opt-
  //     imisation argument, so the corresp "actOpts" bit is set to "true";
  // (*) "OPTIMISE" may be followed by 3 (or none) numerical vals which are
  //     [LoBound InitVal UpBound] of the corresp param; in that case, they
  //     will be used instead of built-in defaults.
  // (*) XXX: For the user convenience, in "Config.ini", the values of "T3"
  //     and "payLoadMass" are given as actual magnitudes (in sec and kg,
  //     resp), whereas the Optimiser always uses relative vals normalised to
  //     [0..1] for the sake of uniformity:
  std::array<bool,NP> actOpts;
  std::vector<double> loBounds;
  std::vector<double> upBounds;
  std::vector<double> initVals;

# ifdef  GetOptParam
# undef  GetOptParam
# endif
# define GetOptParam(I, Type, Name, Scale) \
  Type   Name(NAN);       \
  { \
    std::string paramStr = pt.get<std::string>("Opt." #Name); \
    if (paramStr.substr(0,8) == "OPTIMISE") \
    { \
      /* Parse the params: */ \
      std::vector<std::string> tokens;  \
      boost::split(tokens, paramStr, boost::is_any_of(" \t")); \
      int nTokens = int(tokens.size()); \
      assert(nTokens >= 1); \
      if (nTokens == 2 || nTokens == 3 || nTokens > 4 || \
          tokens[0] != "OPTIMISE")  \
        throw std::invalid_argument \
          ("Ascent21b::FindOptimalAscentCtls: Invalid Params: " + paramStr); \
      /* So yes, the "I"s variable is an optimisation arg: */ \
      actOpts [I] =  true;   \
      if (nTokens == 4)      \
      { \
        /* Set the Range and the InitVal: */ \
        LoBounds[I] = std::atof(tokens[1].data()); \
        InitVals[I] = std::atof(tokens[2].data()); \
        UpBounds[I] = std::atof(tokens[3].data()); \
      } \
      /* For the Dimensioned variable, need to apply the Scale: */ \
      Name  =  InitVals[I] * Scale;    \
      /* The following vectors will be passed to the optimiser: */ \
      loBounds.push_back(LoBounds[I]); \
      upBounds.push_back(UpBounds[I]); \
      initVals.push_back(InitVals[I]); \
    }    \
    else \
    { \
      /* NB: The value for "Name" already has the right magnitude: */ \
      Name       = Type(std::atof(paramStr.data())); \
      actOpts[I] = false; \
    } \
  } \
  assert(IsFinite(Name));

  GetOptParam( 0, double, aAoAHat2,    1.0)
  GetOptParam( 1, double, bAoAHat2,    1.0)
  GetOptParam( 2, Time,   T3,          SC::Stage3BurnDur)
  GetOptParam( 3, double, aAoAHat1,    1.0)
  GetOptParam( 4, double, bAoAHat1,    1.0)
  GetOptParam( 5, Mass,   payLoadMass, massScale)

  // So: How many Optimisation Args have we got?
  int np = 0;
  for (int i = 0; i < NP; ++i)
    np += int(actOpts[size_t(i)]);

  assert(int(loBounds.size()) == np && int(upBounds.size()) == np &&
         int(initVals.size()) == np);

  // Limits for the Constained Variables:
  VelK     maxStartV (pt.get<double>("Opt.MaxStartV" ));
  Pressure QLimit    (pt.get<double>("Opt.QLimit"    ));
  Pressure sepQLimit (pt.get<double>("Opt.SepQLimit" ));
  double   longGLimit(pt.get<double>("Opt.LongGLimit"));
  // Soyuz-2.1b: the limit on the End Time Miss, normally the same as in
  // "FindBurnTime3" ("EndTimeTol"):
  Time     maxEndTimeMiss(pt.get<double>("Opt.MaxEndTimeMiss"));

  proto.ModifyLVParams
  (
    // Mass Params:
    payLoadMass,
    // Ctl Params:
    aAoAHat2,    bAoAHat2,  T3,
    aAoAHat1,    bAoAHat1
  );

  // The status of the result (no optimisation if np=0):
  char const* status = "NOT OPTIMISED";

  if (np != 0)
  {
    //-----------------------------------------------------------------------//
    // Create the NOMAD Optimiser:                                           //
    //-----------------------------------------------------------------------//
    // The Main NOMAD obj:
    NOMAD::MainStep opt;

    // Create and set the NOMAD params:
    // There are 5 Constraints: (StartV, MaxQ, SepQ, MaxLongG, EndTimeMiss),
    // but their actual vals are not required yet:
    constexpr int NC = 5;
    std::shared_ptr<NOMAD::AllParameters> params =
      MkNOMADParams(initVals,    loBounds, upBounds,   NC,
                    optMaxEvals, optSeed,  stopIfFeasible, useVNS, useMT);

    opt.setAllParameters(params);

    std::unique_ptr<NOMADEvaluator> ev
      (new NOMADEvaluator
      (
        params->getEvalParams(),
        &proto,    actOpts,  massScale,
        maxStartV, QLimit,   sepQLimit, longGLimit, maxEndTimeMiss
      ));
    opt.setEvaluator(std::move(ev));

    //-----------------------------------------------------------------------//
    // RUN the Optimiser:                                                    //
    //-----------------------------------------------------------------------//
    opt.start();
    opt.run();
    opt.end();

    //-----------------------------------------------------------------------//
    // Extract the Results:                                                  //
    //-----------------------------------------------------------------------//
    std::vector<NOMAD::EvalPoint> feasPts;
    (void) NOMAD::CacheBase::getInstance()->findBestFeas(feasPts);
    status = "FEASIBLE";

    if (feasPts.empty())
    {
      // Soyuz-2.1b: no result, unless the least infeasible point is requested
      // (for the Final Run and the output only):
      if (!withInfeasible)
        return std::make_pair(std::nullopt, std::nullopt);

      (void) NOMAD::CacheBase::getInstance()->findBestInf(feasPts);
      if (feasPts.empty())
        return std::make_pair(std::nullopt, std::nullopt);
      status = "INFEASIBLE (NOT CLOSED)";
    }

    // If OK: Put the best solution back to "initVals":
    NOMAD::EvalPoint bestF = feasPts[0];

    for (int unsigned j = 0; j < unsigned(np); ++j)
      (initVals)[j] = bestF[j].todouble();

    // The "optMass" found (either the mininimised StartMass, or the maximised
    // PayLoadMass):
    // Mass optMass(bestF.getF(NOMAD::defaultFHComputeType).todouble());
  }
  // Otherwise (if np=0), no optimisation is performed, but the post-processing
  // still is...

  //-------------------------------------------------------------------------//
  // Post-Processing:                                                        //
  //-------------------------------------------------------------------------//
  // In order to get the final "OptRes",  put  the "initVals" (now containing
  // the opt args found) into the "proto" obj, and extract the "optRes"  from
  // the "proto":
  //
  for (int i = 0, j = 0; i < NP; ++i)
  {
    if (!actOpts[size_t(i)])
      continue;
    assert(j < np);

    // XXX: The previously-declared variables are used here. They either hold
    // fixed values from Config.ini, or assigned the optimal vals found:
    switch (i)
    {
    case  0: aAoAHat2    = initVals[unsigned(j)]; break;
    case  1: bAoAHat2    = initVals[unsigned(j)]; break;
    case  2: T3          = initVals[unsigned(j)] * SC::Stage3BurnDur; break;
    case  3: aAoAHat1    = initVals[unsigned(j)]; break;
    case  4: bAoAHat1    = initVals[unsigned(j)]; break;
    case  5: payLoadMass = initVals[unsigned(j)] * massScale;         break;
    default: assert(false);
    }
    ++j;
  }

  proto.ModifyLVParams
  (
    // Mass Params:
    payLoadMass,
    // Ctl Params:
    aAoAHat2,    bAoAHat2,  T3,
    aAoAHat1,    bAoAHat1
  );

  // The Optimisation Result:
  std::optional<OptRes> optRes = OptRes(proto, status);

  if (withFinalRun)
  {
    //-----------------------------------------------------------------------//
    // Perform the Final Run on the optimised params:                        //
    //-----------------------------------------------------------------------//
    // "proto" is already set up with those params, except for the LogLevel:
    //
    proto.m_logLevel    = finalRunLogLevel;
    Base::RunRes runRes = proto.Run();
    return std::make_pair(optRes, runRes);
  }
  // Otherwise, there is no "runRes":
  return std::make_pair(optRes, std::nullopt);
}

//===========================================================================//
// Soyuz-2.1b: "FindOptimalAscentCtls" from the Config.ini File:             //
//===========================================================================//
std::pair<std::optional<Ascent21b::OptRes>,
          std::optional<Ascent21b::Base::RunRes>>
Ascent21b::FindOptimalAscentCtls
(
  std::string const& a_config_ini,
  std::ostream*      a_os       // May be NULL
)
{
  boost::property_tree::ptree  pt;
  boost::property_tree::ini_parser::read_ini(a_config_ini, pt);
  return FindOptimalAscentCtls(pt, a_os);
}
}
// End namespace SpaceBallistics
