// vim:ts=2:et
//===========================================================================//
//                    "Src/Missions/Ascent21b-Integr.cpp":                   //
//    Ascent-to-Orbit of the "Soyuz-2.1b" LV: Bwd Trajectory Integration     //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent21b.h"
#include "LVBase.hpp"
#include "SpaceBallistics/PhysEffects/LVAeroDyn.hpp"
#include "SpaceBallistics/Maths/RKF5.hpp"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage1_Booster.h"

namespace SpaceBallistics
{
namespace
{
  //=========================================================================//
  // Soyuz-2.1b: Local Consts:                                               //
  //=========================================================================//
  // NB: The Stage Models require a defined Launch Instant, but only for the
  // time-stamps of their ECOS vectors; the dynamics does not depend on it, so
  // J2000.0 is used:
  constexpr FlightTime LiftOffFT(0.0_sec, TT(0.0_sec));

  // NB: The Head Unit is the Large (ST) Fairing and the Large InterStage of
  // "Soyuz21b_Head"; the InterStage stays with Stage3 up to the Orbital Ins-
  // ertion:
  constexpr Mass       FairingMass    = SH::FairingLargeMass;
  constexpr Mass       InterStageMass = SH::InterStageLargeMass;
  constexpr Len        FairingD       = SH::FairingLargeMaxD;

  // At the top of a Segment (the later Fwd Time end), the Stage Models  are
  // taken at the left limit, so that the whole Segment sees the same branch:
  constexpr Time       SegEndEps      = Time(1e-9);

  using B = Soyuz21b_Stage1_Booster<'B'>;
}

//===========================================================================//
// "Ascent21b": Non-Default Ctor:                                            //
//===========================================================================//
// Initialises the LV and Mission Params; the Ctl Params are set to the default
// vals:
Ascent21b::Ascent21b
(
  // Stage3:
  Time                    a_T3,
  // Stages 2..3:
  Angle_deg               a_max_aoa2,
  // Stage1:
  Angle_deg               a_max_aoa1,

  // Over-All:
  Mass                    a_payload_mass,

  // Mission Params:
  LenK                    a_h_perigee,
  LenK                    a_h_apogee,
  Angle_deg               a_incl,
  Soyuz2_LaunchPad const& a_pad,

  // Integration/Output Params:
  Time                    a_ode_integr_step,
  std::ostream*           a_os,
  int                     a_log_level
)
: Base
  (
    // Soyuz-2.1b: the Stage1 Params are those of the 4 Boosters of the Stage
    // Model (for info only: the dynamics uses the Stage Models themselves):
    4.0 * B::FullMass,
    double((B::FuelMass + B::OxidMass) / B::FullMass),
    double((B::FuelRem  + B::OxidRem)  / (B::FuelMass + B::OxidMass)),
    B::IspEngSL,
    B::IspEngVac,
    4.0 * To_Len_km(B::ThrustEngVac),
    B::ShutDownThrottlLevel,
    FairingD,
    // Integration and Output Params:
    a_ode_integr_step,
    a_os,
    a_log_level
  ),
  //-------------------------------------------------------------------------//
  // Over-All and some Stage1 Params:                                        //
  //-------------------------------------------------------------------------//
  // Non-Const Over-All Params (may be updated in the course of Optimisation):
  m_payLoadMass   (a_payload_mass),

  //-------------------------------------------------------------------------//
  // Mission Params:                                                         //
  //-------------------------------------------------------------------------//
  m_Rins          (0.0_km), // Initialised below
  m_Vins          (0.0),    //

  //-------------------------------------------------------------------------//
  // Flight Ctl Params:                                                      //
  //-------------------------------------------------------------------------//
  // Stage2:
  m_maxAoA2       (To_Angle(a_max_aoa2)),
  m_aAoA2         (0.0),
  m_bAoA2         (0.0),
  m_aAoAHat2      (0.0),
  m_bAoAHat2      (0.0),

  // Stage1:
  m_maxAoA1       (To_Angle(a_max_aoa1)),
  m_aAoA1         (0.0),
  m_bAoA1         (0.0),
  m_aAoAHat1      (0.0),
  m_bAoAHat1      (0.0),

  //-------------------------------------------------------------------------//
  // Transient Data:                                                         //
  //-------------------------------------------------------------------------//
  m_mode          (FlightMode::Burn3),
  m_fairingSepTime(NAN),
  m_eventStr      (),
  m_nextOutputTime(),
  m_maxQ          (0.0),
  m_sepQ          (0.0),
  m_maxLongG      (0.0),

  //-------------------------------------------------------------------------//
  // Soyuz-2.1b: the Stage3 BurnTime, the Segments, the Stage Models:        //
  //-------------------------------------------------------------------------//
  m_T3            (NAN),    // Set by "SetBurnTime3" below
  m_Tflight       (NAN),    //
  m_T2            (NAN),    //
  m_segTop        (NAN),
  m_segBot        (NAN),
  m_segEvent      (nullptr),
  m_stack         { Time(NAN), Pressure(NAN), Density(NAN), Mass(NAN),
                    MassRate(NAN), ForceK(NAN) }
{
  //-------------------------------------------------------------------------//
  // Checks:                                                                 //
  //-------------------------------------------------------------------------//
  if (!(IsPos(m_maxAoA2) && IsPos(m_maxAoA1)))
    throw std::invalid_argument("Ascent21b::Ctor: Invalid LV Param(s)");

  // XXX: We assume that launch is NOT from the North or South Pole:
  Angle launchLat = a_pad.Latitude();
  if (!(a_h_perigee       >=  100.0_km && a_h_perigee <= a_h_apogee &&
        a_incl            >= 0.0_deg   && a_incl      <= 180.0_deg  &&
        Abs(launchLat)     < PI_2      && !IsNeg(m_payLoadMass)))
    throw std::invalid_argument("Ascent21b::Ctor: Invalid Mission Param(s)");

  //-------------------------------------------------------------------------//
  // Mission Params:                                                         //
  //-------------------------------------------------------------------------//
  // Semi-Major Axis and the Perigee radius-vector:
  LenK a = R + 0.5 * (a_h_perigee + a_h_apogee);
  m_Rins = R +        a_h_perigee;

  // The required orbital velocity at the orbital insertion point (which is
  // always assumed to be the Perigee), using the Energy Integral:
  VelK Vq = SqRt(K * (2.0 / m_Rins - 1.0 / a));

  // The Orbit Inclination and the Launch Site Latitude:
  double cosI   = Cos(To_Angle_rad(a_incl));
  double cosLat = Cos(launchLat);
  assert(cosLat > 0.0);

  if (std::fabs(cosI) > cosLat)
    throw std::invalid_argument("Ascent21b::Ctor: Inclination Unreachable");

  // "A" is the Launch Azimuth:
  double sinA   = cosI / cosLat;
  assert(std::fabs(sinA) <= 1.0);

  // We will assume -Pi/2 <= A <= Pi/2, so cosA >= 0:
  double cosA   = SqRt(1.0 - Sqr(sinA));

  // The Earth Rotation Velocity at the Launch Latitude:
  VelK   ERV    =
      TwoPi<double> * R / BodyData<Body::Earth>::SiderealRotationPeriod *
      cosLat;

  // The components of the velocity (in the Launch COS) which must be attained
  // by the LV at the orbital insertion point, taking the Earth rotation into
  // account:
  VelK   Vlat   = Vq * sinA - ERV;
  VelK   Vlong  = Vq * cosA;
  // So the LV velocity at the orbital insertion point must be:
  m_Vins        = SqRt(Sqr(Vlat) + Sqr(Vlong));

  //-------------------------------------------------------------------------//
  // Soyuz-2.1b: the AoA Domains and the Stage3 BurnTime:                    //
  //-------------------------------------------------------------------------//
  // NB: The Stage1 AoA law spans the whole "Burn12" mode,  ie up to the Stage1
  // Separation (not the Cut-Off), see "AoA":
  Base::m_T1 = Time(SC::Stage1SepTime) - Time(SC::LiftOffTime);

  // The Stage3 BurnTime and the Ctls depending on it:
  SetBurnTime3(a_T3);
}

//===========================================================================//
// "Run": Integrate the Ascent Trajectory:                                   //
//===========================================================================//
Ascent21b::Base::RunRes Ascent21b::Run()
{
  //-------------------------------------------------------------------------//
  // Run the integration BACKWARDS from the Orbital Insertion point          //
  //-------------------------------------------------------------------------//
  // (@ t=0):
  // Angular LV velocity at the orbital insertion point:
  AngVel omega0 = 1.0_rad * m_Vins / m_Rins;

  // The initial Radial Velocity is 0, and the final delta of Spent Mass
  // is also 0, and the initial Polar Angle is 0:
  Base::StateV s0 =
    std::make_tuple(m_Rins, VelK(0.0), omega0, 0.0_kg, 0.0_rad);

  // The RHS and the Call-Back Lambdas (IsAscent=true):
  auto rhs =
    [this](Base::StateV const&  a_s, Time a_t) -> Base::DStateV
    { return this->Base::ODERHS(a_s, a_t, true); };

  auto cb  =
    [this]
    (
      Base::StateV*        a_s,
      Time                 a_t,
      Base::DStateV const& a_ds,
      Base::StateV  const& a_prev_s,
      Time                 a_prev_t,
      Base::DStateV const& a_prev_ds
    )
    -> bool
    { return this->ODECB (a_s, a_t, a_ds, a_prev_s, a_prev_t, a_prev_ds); };

  // NB: Transient fields have been initialised in the Ctor. XXX: IMPORTANT:
  // It is currently assumed that "Run" is invoked only once per the object
  // lifetime!
  if (m_mode != FlightMode::Burn3)
    throw std::runtime_error
          ("Ascent21b::Run Repeated invocations not allowed");

  //-------------------------------------------------------------------------//
  // Soyuz-2.1b: the Segments:                                               //
  //-------------------------------------------------------------------------//
  // Instead of the max ascent duration: the Segments delimited by the Stage
  // Models' events (in the decreasing Fwd Time), so that the Stage Models
  // have no discontinuities inside a Segment;  the Fairing Separation is the
  // only event inside a Segment (it is detected by the Call-Back):
  struct Event
  {
    Time        m_tau;
    char const* m_name;
  };
  Event const events[]
  {
    { Time(SC::Stage3AftJetTime), "Stage3 Aft Section Jettison"     },
    { Time(SC::Stage2SepTime),    "Stage2 Separation"               },
    { Time(SC::Stage2CutOffTime), "Stage2 Cut-Off, Stage3 Ignition" },
    { Time(S2::MainCutOffTime),   "Stage2 Main Engine Cut-Off"      },
    { Time(S2::ThrottlTime),      "Stage2 Throttle-Down"            },
    { Time(SC::Stage1SepTime),    "Stage1 Separation"               },
    { Time(SC::Stage1CutOffTime), "Stage1 Cut-Off"                  },
    { Time(B ::ThrottlTime),      "Stage1 Throttle-Down"            },
    { Time(B ::FullThrustTime),   "Stage1 Full Thrust"              },
    { Time(SC::LiftOffTime),      nullptr                           }
  };
  static_assert(SC::LiftOffTime      < B ::FullThrustTime   &&
                B ::FullThrustTime   < B ::ThrottlTime      &&
                B ::ThrottlTime      < SC::Stage1CutOffTime &&
                SC::Stage1CutOffTime < SC::Stage1SepTime    &&
                SC::Stage1SepTime    < S2::ThrottlTime      &&
                S2::ThrottlTime      < S2::MainCutOffTime   &&
                S2::MainCutOffTime   < SC::Stage2CutOffTime &&
                SC::Stage2CutOffTime < SC::Stage2SepTime    &&
                SC::Stage2SepTime    < SC::Stage3AftJetTime);

  // The FlightMode below the Fwd Time "a_tau":
  auto modeBelow =
    [](Time a_tau) -> FlightMode
    {
      return
        (a_tau > Time(SC::Stage2SepTime))
        ? FlightMode::Burn3  :
        (a_tau > Time(SC::Stage2CutOffTime))
        ? FlightMode::Burn23 :
        (a_tau > Time(SC::Stage1SepTime))
        ? FlightMode::Burn2
        : FlightMode::Burn12;
    };

  Time           tEnd = 0.0_sec;      // Set by the Segments below
  m_nextOutputTime    = 0.0_sec;      // Just to make sure

  Base::RunRes   res;                 // "Error" as yet
  try
  {
    //-----------------------------------------------------------------------//
    // Actually Run the Integrator, Segment by Segment!                      //
    //-----------------------------------------------------------------------//
    // XXX: "RKF5" does not re-compute the target time after a step reduction,
    // so the step control is disabled (eps=0): with all Stage Model discont-
    // inuities at the Segment boundaries, the fixed-step RKF5 is accurate
    // enough:
    //
    Time        top   = m_Tflight;
    char const* topEv = nullptr;   // The Orbital Insertion is logged anyway

    for (Event const& ev: events)
    {
      if (ev.m_tau >= top)
        continue;                  // Above the Orbital Insertion

      m_segTop   = top;
      m_segBot   = ev.m_tau;
      m_segEvent = topEv;
      m_mode     = modeBelow(top);

      Time tTop  = top       - m_Tflight;
      Time tBot  = ev.m_tau  - m_Tflight;
      tEnd       =
        RKF5(&s0,  tTop, tBot, rhs,
             -Base::m_odeIntegrStep, -Base::m_odeIntegrStep, 0.0, &cb,
              Base::m_os);
      assert(tEnd >= tBot);

      // Stopped by the Call-Back (H <= 0 or the Lift-Off)?
      if (tEnd != tBot || m_mode == FlightMode::UNDEFINED ||
          !IsPos(std::get<0>(s0) - R))
        break;

      top   = ev.m_tau;
      topEv = ev.m_name;
    }

    // If we got here: Integration has run to completion, but not to the Singu-
    // lar Point. Invoke the "PostProcessRun" (IsAscent=true):
    res  = Base::PostProcessRun(s0, tEnd, true);
  }
  catch (Base::NearSingularityExn const& ns)
  {
    // We have reached a vicinity of the Singular Point:
    // This is a NORMAL (and moreover, a desirable) outcome. Compute a more
    // precise singular point position;  IsAscent=true:
    res  = Base::LocateSingularPoint(ns, true);
  }
  catch (FallingBackExn const& fb)
  {
    // No point in continuing: we are falling back to Earth instead of ascend-
    // ing to orbit:
    VelK V = SqRt(Sqr(fb.m_Vr) + Sqr(fb.m_Vhor));

    if (Base::m_os != nullptr)
#     pragma omp critical(Output)
      (*m_os)
        << "# ERROR: Falling Back to Earth: t="           << fb.m_t.Magnitude()
        << " sec,  h="         << (fb.m_r - R).Magnitude() << " km, l="
        << (R * fb.m_phi / 1.0_rad).Magnitude()            << " km, Vr="
        << fb.m_Vr.Magnitude() << " m/sec, V="             << V.Magnitude()
        << " km/sec, mode="    << ToString(m_mode)         << std::endl;

    // XXX: We do not propagate vals carried by "fb" to "RunRes" which remains
    // an error...
  }
  // XXX: Any other exceptions are propagated to the top level, it's better not
  // to hide them...
  // The result:
  if (Base::m_os != nullptr && Base::m_logLevel >= 3)
#   pragma omp critical(Output)
    *Base::m_os
      << "# t=" << res.m_T.Magnitude() << " sec: Integration Done, RC="
      << Base::ToString(res.m_rc)      << ", Tflight="
      << m_Tflight.Magnitude()         << " sec, T3="
      << m_T3.Magnitude()              << " sec" << std::endl;
  return res;
}

//===========================================================================//
// ODE CallBack (invoked after the completion of each RKF5 step):            //
//===========================================================================//
// HERE FlightMode switching occurs, so this method is non-"const";
// some args are unused yet:
//
bool Ascent21b::ODECB
(
  Base::StateV*   a_s, Time a_t,      Base::DStateV const&,
  Base::StateV const&, Time a_prev_t, Base::DStateV const&
)
{
  assert(a_s != nullptr && !IsPos(a_t));
  LenK   r             = std::get<0>(*a_s);
  LenK   h             = r - R;             // Altitude
  VelK   Vr            = std::get<1>(*a_s);
  AngVel omega         = std::get<2>(*a_s);
  VelK   Vhor          = r * omega / 1.0_rad;
  Mass   spentPropMass = std::get<3>(*a_s);
  Angle  phi           = std::get<4>(*a_s);

  // NB: "spentPropMass" decreases over time, so increases in Bwd time.
  // It is 0 at Orbit Insertion Time (t=0):
  assert(IsPos(r) && !IsNeg(spentPropMass));

  //-------------------------------------------------------------------------//
  // Check if we are approaching the Singularity:                            //
  //-------------------------------------------------------------------------//
  // XXX: Unlike Base::ODERHS, here we set the angular velocity to 0 if "Vhor"
  // (not |Vhor|!) is below the threshold, because we KNOW that this may only
  // happen near the Singular Point:
  //
  if (Vhor < Base::SingVhor)
  {
    omega              = AngVel();
    std::get<2>(*a_s)  = AngVel();
    Vhor               = VelK  ();
  }
  if (IsZero(omega) && Abs(Vr) < Base::SingV)
    // This is not an error -- just stop integration now. XXX: Small negative
    // vals of "Vr" are handled in this case as well. In theory, V=0 does not
    // always imply a singularity; but in practice, it does:
    //
    throw Base::NearSingularityExn(r, Vr, spentPropMass, phi, a_t);

  // Otherwise, we allow Vr < 0 ONLY during "Burn3" and if the Ascent traject-
  // ory is above the orbital insertion point; in all other cases, it means
  // that we are falling back to Earth, so stop integration immediately:
  if (IsNeg(Vr) && !(m_mode == FlightMode::Burn3 && r > m_Rins))
  {
    // If it is small (< 10 m/sec), just correct it:
    if (Abs(Vr) < VelK(0.01))
      Vr = VelK(0.0);
    else
      // Otherwise, throw an exception; when handled it will invalidate this
      // solution:
      throw FallingBackExn(r, Vr, Vhor, spentPropMass, phi, a_t);
  }

  //-------------------------------------------------------------------------//
  // Generic Case:                                                           //
  //-------------------------------------------------------------------------//
  // Compute the Non-Gravitational Acceleration Components:
  Mass          m      = LVMass(*a_s, a_t);  // Using the "old" Mode yet!
  VelK          V;
  Angle         psi, aoa;
  EAM::AtmConds atm;
  MassRate      burnRate;
  ForceK        thrust;
  double        lvAxis[2];
  AccK          ngAcc [2];
  double        longG;

  Base::NonGravForces
    (a_t, r, Vr, Vhor, m, &V, &psi, &aoa, &atm, &burnRate, &thrust,
     lvAxis, ngAcc, &longG);

  LenK     L   = R * phi / 1.0_rad; // Down-Range Earth Distance, <= 0
  // AeroDynamic Conditions:
  double   M   = Base::Mach(atm, V);
  Pressure pa  = std::get<0>(atm);
  Density  rho = std::get<1>(atm);
  auto     V2  = To_Len_m(Sqr(V));
  Pressure Q   = 0.5 * rho * V2;

  // Curr event for logging (if any):
  m_eventStr.clear();

  if (IsZero(a_t))
    m_eventStr = "Stage3 Cut-Off, Orbital Insertion";

  //-------------------------------------------------------------------------//
  // Soyuz-2.1b: the Stage Models' Events and the Lift-Off:                  //
  //-------------------------------------------------------------------------//
  // The modes Burn3 .. Burn12 are switched by "Run" at the Segment bounds; an
  // event is logged at the top of the curr Segment:
  if (m_segEvent != nullptr && a_t == m_segTop - m_Tflight)
  {
    m_eventStr = m_segEvent;

    // XXX: Stage1 Separation: record the corresp Q:
    if (m_segTop == Time(SC::Stage1SepTime))
      m_sepQ   = Q;
  }

  // Switching to UNDEFINED @ the Lift-Off Time:
  if (a_t <= - m_Tflight)
  {
    m_mode     = FlightMode::UNDEFINED;
    m_eventStr = "Lift-Off";
  }

  //-------------------------------------------------------------------------//
  // In any mode, detect the Fairing Separation Condition:                   //
  //-------------------------------------------------------------------------//
  // In the reverse time, it's when the Dynamic Pressure becomes HIGHER
  // that the threshold; Soyuz-2.1b: but the instant is kept within the range
  // given by "Soyuz21b_Consts":
  Time tau = a_t + m_Tflight;   // Fwd Flight Time

  if (!IsFinite(m_fairingSepTime) && tau <= Time(SC::FairingJetTimeMax) &&
      (Q >= FairingSepCond        || tau <= Time(SC::FairingJetTimeMin)))
  {
    m_fairingSepTime = a_t;
    m_eventStr       = "Fairing Separation";
  }

  //-------------------------------------------------------------------------//
  // Monitor the Constraints:                                                //
  //-------------------------------------------------------------------------//
  // If there was any event, re-compute the LV Mass, Forces and Accelerations:
  if (!m_eventStr.empty())
  {
    m    = LVMass(*a_s, a_t);
    Base::NonGravForces
      (a_t, r, Vr, Vhor, m, &V, &psi, &aoa, &atm, &burnRate, &thrust,
      lvAxis, ngAcc, &longG);
  }
  // Dynamic Pressure:
  m_maxQ     = std::max(Q, m_maxQ);

  // LongG:
  assert(!IsNeg(longG));
  m_maxLongG = std::max(longG,  m_maxLongG);

  //-------------------------------------------------------------------------//
  // Integration Stopping Conds:                                             //
  //-------------------------------------------------------------------------//
  bool cont = true;

  // Do NOT continue if:
  // (*) we are on the surface (stopping by "h");
  // (*) we are in the UNDEFINED mode (stopping by the Lift-Off Time);
  // NB: both are NORMAL termination conds (as opposed to throwing
  //     the "StopNowExn"):
  if (!IsPos(h))
  {
    cont        = false;
    m_eventStr += ": Integration Stopped @ H <= 0";
  }
  if (m_mode == FlightMode::UNDEFINED)
  {
    cont        = false;
    m_eventStr += ": Integration Stopped @ Lift-Off Time";
  }

  //-------------------------------------------------------------------------//
  // Log the Curr Event (if any):                                            //
  //-------------------------------------------------------------------------//
  // Angles in Degrees:
  Angle_deg psi_deg = To_Angle_deg(psi);
  Angle_deg aoa_deg = To_Angle_deg(aoa);

  if (!m_eventStr.empty() && Base::m_os != nullptr && Base::m_logLevel >= 3)
#   pragma omp critical(Output)
    *Base::m_os
      << "# t="           << a_t.Magnitude()
      << " sec: "         << m_eventStr
      << ": Mode="        << ToString(m_mode)
      << ", h="           << h.Magnitude()
      << " km, L="        << L.Magnitude()
      << " km, V="        << V.Magnitude()
      << " km/sec, Vr="   << Vr.Magnitude()
      << " km/sec, Vhor=" << Vhor.Magnitude()
      << " km/sec, psi="  << psi_deg.Magnitude()
      << " deg, AoA="     << aoa_deg.Magnitude()
      << " deg, LVMass="  << m.Magnitude()
      << " kg, Q="        << Q.Magnitude()
      << " Pa, M="        << M
      << ", LongG="       << longG
      << std::endl;

  //-------------------------------------------------------------------------//
  // Main Output:                                                            //
  //-------------------------------------------------------------------------//
  // Occurs with a 0.1 msec step, or if we are going to stop now;
  // Soyuz-2.1b: plus the 18th column, the Spent Propellant Mass (kg):
  if (Base::m_os != nullptr && Base::m_logLevel >= 4 &&
     (!cont || a_t <= m_nextOutputTime))
  {
    m_nextOutputTime -= 0.1_sec;   // NB: Bwd Integration!

    // Thrust is more conveniently reported in kgf:
    auto tkg = thrust / g0K;
    Time dt  = a_t - a_prev_t;
    assert(!IsPos(dt));            // NB: dt=0 at the start of each Segment

#   pragma omp critical(Output)
    *Base::m_os
      << a_t.Magnitude()     << '\t' << h.Magnitude()        << '\t'
      << L.Magnitude()       << '\t' << Vr.Magnitude()       << '\t'
      << Vhor.Magnitude()    << '\t' << V.Magnitude()        << '\t'
      << psi_deg.Magnitude() << '\t' << aoa_deg.Magnitude()  << '\t'
      << m.Magnitude()       << '\t' << ToString(m_mode)     << '\t'
      << tkg.Magnitude()     << '\t' << burnRate.Magnitude() << '\t'
      << Q.Magnitude()       << '\t' << M                    << '\t'
      << longG               << '\t' << pa.Magnitude()       << '\t'
      << dt.Magnitude()      << '\t' << spentPropMass.Magnitude()
      << std::endl;
  }
  return cont;
}

//===========================================================================//
// Angle-of-Attack and Thrust Vector Elevation:                              //
//===========================================================================//
//        AoA   theta
std::pair<Angle,Angle> Ascent21b::AoA(Time a_t, Angle a_psi) const
{
  // XXX: We do not check for sin(a_psi) < 0 here (ie switching from Ascent to
  // Descent motion),  as this would require a wider context than available in
  // this method. Such checks are performed in "ODECB".
  //
  // Also, the "retrograde" movement (|psi| > Pi/2, and thus omega < 0,  since
  // omega ~ Vhor ~ cos(psi)) may in general occur due to rounding errors or
  // integration instabilities. It is difficult to say precisely which values
  // of "psi" are accepatble and which are not, so we just log such events:
  //
  if (Abs(a_psi) > 1.001 * PI_2 &&
      Base::m_os != nullptr     && Base::m_logLevel >= 1)
#   pragma omp critical(Output)
    *Base::m_os << "# Ascent21b::AoA: WARNING: psi="
                << To_Angle_deg(a_psi).Magnitude() << " deg" << std::endl;

  Angle  aoa = 0.0_rad;

  switch (m_mode)
  {
  case FlightMode::Burn12:
  {
    // Soyuz-2.1b: the Stage1 law: "tau" is the time since the Lift-Off,  so
    // AoA1(tau=0)=0, ie @ launch; 0 <= tau <= T1:
    Time tau = std::min(std::max(a_t + m_Tflight, 0.0_sec), Base::m_T1);
    aoa      = std::max(tau * (m_aAoA1 * tau + m_bAoA1), 0.0_rad);
    if (aoa > 1.01 * m_maxAoA1)
      throw std::logic_error
            ("Burn12: MaxAoA1=" + std::to_string(m_maxAoA1.Magnitude()) +
             " exceeded: AoA="  + std::to_string(aoa      .Magnitude()));
    break;
  }

  case FlightMode::Burn2:
  case FlightMode::Burn23:
  case FlightMode::Burn3:
  {
    // Arg: "t" is time from Stage3 cut-off (same as orbital insertion time),
    // so it is just the main time "a_t": -T2 <= t <= 0; enforce those constr-
    // aints. NB: Formulas similar to those for Stage2, but with the INVERTED
    // sign of "b"; AoA2(t=0)=0, ie @ orbital insertion:
    assert(!IsPos(a_t));
    Time t = std::max(a_t, -m_T2);
    aoa    = std::max(t *  (m_aAoA2 * t - m_bAoA2), 0.0_rad);
    if (aoa > 1.01 * m_maxAoA2)
      throw std::logic_error
            ("Burn2..3: MaxAoA2=" + std::to_string(m_maxAoA2.Magnitude()) +
             " exceeded: AoA="    + std::to_string(aoa      .Magnitude()));
    break;
  }
  default: ;
  }
  // XXX: For the Ascent, we can always assume 0 <= aoa <= Pi/2:
  assert(!IsNeg(aoa) && aoa <= PI_2);

  // XXX: In this case, we assume there is no chamber/nozzle gimbaling,  and
  // therefore, theta = psi + AoA; but "theta" is limited by Pi/2, otherwise
  // the Thrust vector will point towards negative polar angles. Yet theta < 0
  // is allowed under some circumstances (because Vr < 0 may be allowed):
  Angle theta = a_psi + aoa;

  // XXX: We apply the following limits on "theta" only, since "aoa" is already
  // in the valid range:
  theta = std::max(-PI_2, std::min(theta, PI_2));
  return  std::make_pair(aoa, theta);
}

//===========================================================================//
// Atmospheric Conditions and Aerodynamic Drag and Lift Forces:              //
//===========================================================================//
//         Conds          Drag      Lift
std::tuple<EAM::AtmConds, ForceK,   ForceK>
Ascent21b::AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA) const
{
  auto atm  = EAM::GetAtmConds(std::max(a_r - R, 0.0_km));
  double M  = Base::Mach(atm, a_v);
  assert(IsFinite(M) &&  M >= 0.0);

  // Soyuz-2.1b: rather than one "Base::m_crosS", the reference area is that
  // of the curr configuration: the disc of the widest body on the axis (Fair-
  // ing, Stage2 or Stage3)  plus those of the 4 Boosters, overlaps neglected.
  // The Fairing is attached iff its Separation has been passed in the Bwd
  // time:
  bool  boosters = m_mode == FlightMode::Burn12 ||
                   m_mode == FlightMode::UNDEFINED;
  bool  core     = boosters || m_mode == FlightMode::Burn2 ||
                               m_mode == FlightMode::Burn23;
  bool  fairing  = IsFinite(m_fairingSepTime);
  Len   axisD    =
    fairing ? std::max(FairingD, S2::MaxD) : core ? S2::MaxD : S3::D;
  AreaK crosS    = 0.25 * Pi<double> * Sqr(To_Len_km(axisD));
  if (boosters)
    crosS       += Pi<double> * Sqr(To_Len_km(B::MaxD));

  // The Aerodynamic Force Main Term:
  ForceK F  = 0.5 * To_Len_km(std::get<1>(atm)) * Sqr(a_v) * crosS;

  // The Drag and Lift Coeffs, using the default model:
  double cD = LVAeroDyn::cD(M, a_AoA);
  double cL = LVAeroDyn::cL(M, a_AoA);

  return std::make_tuple(atm, cD * F, cL * F);
}

//===========================================================================//
// Soyuz-2.1b: the Vehicle from the Stage Models:                            //
//===========================================================================//
// Rather than the "model" LV of "Ascent2", the LV Mass, the Propellant Burn-
// Rate and the Thrust are the sum of the Stage Models ("GetStack") at the Fwd
// Flight Time within the curr Segment ("Tau"), at the actual ambient pressure:
//
//===========================================================================//
// "Tau": The Fwd Flight Time within the curr Segment:                       //
//===========================================================================//
// NB: Clamped to [m_segBot, m_segTop - SegEndEps]:  at the bottom, rounding of
// (t + Tflight) must not switch the Stage Models to the branch below; and the
// Singular Point may be located slightly below the bottom (the Lift-Off):
Time Ascent21b::Tau(Time a_t) const
{
  Time tau = a_t + m_Tflight;
  return std::max(std::min(tau, m_segTop - SegEndEps), m_segBot);
}

//===========================================================================//
// "GetStack": The Stage Models' Sum:                                        //
//===========================================================================//
Ascent21b::Stack const&
Ascent21b::GetStack(Time a_tau, Pressure a_p, Density a_rho) const
{
  assert(!(IsNeg(a_tau) || IsNeg(a_p) || IsNeg(a_rho)));

  // NB: "Soyuz21b_Stage2" asserts p <= P0, so guard against rounding @ h=0:
  Pressure p = std::min(a_p, EAM::P0);

  if (a_tau == m_stack.m_tau && p == m_stack.m_p)
    return m_stack;

  // A Stage contributes its mass and thrust until its separation (the Stage
  // Models themselves give 0 thrust after Cut-Off):
  bool       boosters = a_tau < Time(SC::Stage1SepTime);
  bool       core     = a_tau < Time(SC::Stage2SepTime);
  FlightTime ft       = LiftOffFT + a_tau;

  Mass       mass     = m_payLoadMass + InterStageMass;
  MassRate   massDot  (0.0);
  Force      thrust   (0.0);

  auto add =
    [&mass, &massDot, &thrust](StageDynParams<LVSC::Soyuz21b> const& a_dp)
    -> void
    {
      mass    += a_dp.m_fullMass;
      massDot += a_dp.m_fullMassDot;
      thrust  += a_dp.m_thrust[0];
    };

  // XXX: The Booster's AeroFin model is given ZERO velocity and wind (it
  // asserts against a 90 deg plate AoA otherwise); the AeroFin force is not
  // used in the "StageDynParams" anyway:
  if (boosters)
  {
    constexpr B::VernDeflections vd0;
    constexpr ME::VelVE          v0;
    constexpr Angle_deg          fin0;
    add(Soyuz21b_BlockB::GetDynParams(ft, vd0, p, a_rho, v0, v0, fin0));
    add(Soyuz21b_BlockV::GetDynParams(ft, vd0, p, a_rho, v0, v0, fin0));
    add(Soyuz21b_BlockG::GetDynParams(ft, vd0, p, a_rho, v0, v0, fin0));
    add(Soyuz21b_BlockD::GetDynParams(ft, vd0, p, a_rho, v0, v0, fin0));
  }
  if (core)
  {
    constexpr S2::VernDeflections vd0;
    add(S2::GetDynParams(ft, vd0, p));
  }
  S3::ChamberDeflections const cd0 {};
  add(S3::GetDynParams(ft, cd0));

  assert(IsPos(mass) && !IsPos(massDot) && !IsNeg(thrust));
  m_stack = Stack
  {
    .m_tau      = a_tau,
    .m_p        = p,
    .m_rho      = a_rho,
    .m_mass     = mass,
    .m_burnRate = - massDot,
    .m_thrust   = To_Len_km(thrust)
  };
  return m_stack;
}

//===========================================================================//
// Propellant Burn Rate: May be variable:                                    //
//===========================================================================//
MassRate Ascent21b::PropBurnRate
(
  Mass UNUSED_PARAM(a_m),
  LenK a_h,
  Time a_t
)
const
{
  auto atm = EAM::GetAtmConds(std::max(a_h, 0.0_km));
  return GetStack(Tau(a_t), std::get<0>(atm), std::get<1>(atm)).m_burnRate;
}

//===========================================================================//
// Thrust:                                                                   //
//===========================================================================//
// NB: "LVBase" invokes "Thrust" right after "PropBurnRate" for the same ins-
// tant, so the Fwd Time of the latter (cached in "m_stack") is used:
//
ForceK Ascent21b::Thrust(MassRate DEBUG_ONLY(a_burn_rate), Pressure a_p)
const
{
  assert(!(IsNeg(a_burn_rate) || IsNeg(a_p)));
  Stack const& st = GetStack(m_stack.m_tau, a_p, m_stack.m_rho);
  assert(IsZero(a_burn_rate) == IsZero(st.m_burnRate));
  return st.m_thrust;
}

//===========================================================================//
// Current LV Mass (incl the PayLoad and Fairing if attached):               //
//===========================================================================//
// NB: From the Stage Models, not from the SpentPropMass (which is integrated
// independently, and must be consistent with the former):
//
Mass Ascent21b::LVMass(Base::StateV const& a_s, Time a_t) const
{
  auto atm = EAM::GetAtmConds(std::max(std::get<0>(a_s) - R, 0.0_km));
  Mass m   = GetStack(Tau(a_t), std::get<0>(atm), std::get<1>(atm)).m_mass;

  // Add the FairingMass if the Fairing has NOT separated @  "a_t". In the
  // Bwd time, it means that "m_fairingSepTime" is known and "a_t" is below
  // it:
  if (IsFinite(m_fairingSepTime) && a_t < m_fairingSepTime)
    m += FairingMass;

  // All Done:
  assert(IsPos(m));
  return m;
}

//===========================================================================//
// Soyuz-2.1b: the Ctl Params and the Stage3 BurnTime:                       //
//===========================================================================//
// NB: Unlike "Ascent2" ("Ascent2-Opt.cpp"), "SetCtlParams" is here, as it is
// also used without the Optimiser (by "FindBurnTime3"); the BurnRates are not
// Ctls, and the AoA coeffs depend on "T2", ie on the Stage3 BurnTime:
//
//===========================================================================//
// "SetCtlParams":                                                           //
//===========================================================================//
void Ascent21b::SetCtlParams
(
  double   a_aAoAHat2,  double   a_bAoAHat2,
  double   a_aAoAHat1,  double   a_bAoAHat1
)
{
  //-------------------------------------------------------------------------//
  // Check what we got:                                                      //
  //-------------------------------------------------------------------------//
  if (!(0.0 <= a_aAoAHat2    && a_aAoAHat2 <= 1.0  &&
        0.0 <= a_bAoAHat2    && a_bAoAHat2 <= 1.0  &&
        0.0 <= a_aAoAHat1    && a_aAoAHat1 <= 1.0  &&
        0.0 <= a_bAoAHat1    && a_bAoAHat1 <= 1.0))
      throw std::invalid_argument("Ascent21b::SetCtlParams: Invalid Param(s)");

  m_aAoAHat2     = a_aAoAHat2;
  m_bAoAHat2     = a_bAoAHat2;
  m_aAoAHat1     = a_aAoAHat1;
  m_bAoAHat1     = a_bAoAHat1;
  SetAoACoeffs();
}

//===========================================================================//
// "SetAoACoeffs":                                                           //
//===========================================================================//
void Ascent21b::SetAoACoeffs()
{
  //-------------------------------------------------------------------------//
  // AoA Coeffs:                                                             //
  //-------------------------------------------------------------------------//
  // Stage2:
  m_bAoA2        = 4.0 * m_maxAoA2 / m_T2 * m_bAoAHat2;
  AngAcc aAoALo2 = - m_bAoA2       / m_T2;
  AngAcc aAoAUp2 =
    (m_bAoAHat2  < 0.5)
    ? (m_maxAoA2     / m_T2 - m_bAoA2) / m_T2
    : - Sqr(m_bAoA2) / (4.0 * m_maxAoA2);
  m_aAoA2        =     (1.0 - m_aAoAHat2) * aAoALo2 + m_aAoAHat2 * aAoAUp2;

  // Stage1:
  m_bAoA1        = 4.0 * m_maxAoA1 / Base::m_T1 * m_bAoAHat1;
  AngAcc aAoALo1 = - m_bAoA1       / Base::m_T1;
  AngAcc aAoAUp1 =
    (m_bAoAHat1  < 0.5)
    ? (m_maxAoA1     / Base::m_T1  - m_bAoA1) / Base::m_T1
    : - Sqr(m_bAoA1) / (4.0 * m_maxAoA1);
  m_aAoA1        =     (1.0 - m_aAoAHat1) * aAoALo1 + m_aAoAHat1 * aAoAUp1;
}

//===========================================================================//
// "SetBurnTime3":                                                           //
//===========================================================================//
void Ascent21b::SetBurnTime3(Time a_T3)
{
  // NB: Stage3 must burn alone at the Orbital Insertion (ie beyond the Stage2
  // Separation), and not longer than the Stage Model allows:
  if (!(a_T3 >  Time(SC::Stage2SepTime) - Time(SC::Stage2CutOffTime) &&
        a_T3 <= SC::Stage3BurnDur))
    throw std::invalid_argument("Ascent21b::SetBurnTime3: Invalid T3");

  m_T3      = a_T3;
  m_Tflight = Time(SC::Stage2CutOffTime) - Time(SC::LiftOffTime) + m_T3;
  m_T2      = m_Tflight - Base::m_T1;
  m_segTop  = m_Tflight;
  m_segBot  = 0.0_sec;
  assert(IsPos(m_T2));

  // The AoA2 coeffs depend on "T2":
  SetAoACoeffs();
}

//===========================================================================//
// "GetReserve3":                                                            //
//===========================================================================//
Mass Ascent21b::GetReserve3() const
{
  S3::ChamberDeflections const cd0 {};
  auto dp = S3::GetDynParams
            (LiftOffFT + (Time(SC::Stage3IgnTime) - Time(SC::LiftOffTime) +
                          m_T3), cd0);
  return (dp.m_fuelMass - S3::FuelRem) + (dp.m_oxidMass - S3::OxidRem);
}
}
// End namespace SpaceBallistics
