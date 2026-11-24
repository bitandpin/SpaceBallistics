// vim:ts=2:et
//===========================================================================//
//                   "SpaceBallistics/Missions/Ascent21b.h":                 //
//        Ascent-to-Orbit of the "Soyuz-2.1b" LV, using its Stage Models     //
//===========================================================================//
// "Ascent21b" is the method of "Ascent2" (BACKWARD integration from the Orb-
// ital Insertion at t=0 to the Lift-Off, via "LVBase") applied to the Stage
// Models of the "Soyuz-2.1b" LV; read "Ascent2" first. These files follow the
// "Ascent2" ones section by section; the differences are gathered in blocks
// marked "Soyuz-2.1b":
// (*) the Vehicle from the Stage Models: the LV Mass, BurnRate and Thrust are
//     the sum of the Stage Models at the Fwd Flight Time tau = t + Tflight and
//     the actual ambient pressure, rather than a "model" LV, so the BurnRates
//     are not Ctl Params;
// (*) 3 Stages with parallel staging: the FlightModes "Burn12" .. "Burn3" and
//     the integration by Segments between the Stage Models' events (Lift-Off,
//     Separations, ...); no Ballistic Gap and no Stage1 return;
// (*) the Stage3 BurnTime closure: "T3" (the Block I is cut off by guidance
//     when the orbit is reached) is the closing param, such that  the Bwd
//     integration ends at the Lift-Off time; "FindBurnTime3" searches for it,
//     and the Optimiser treats it as an Opt Param with an extra Constraint.
//
#pragma  once
#include "SpaceBallistics/Missions/LVBase.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/LaunchPads.hpp"
#include <boost/property_tree/ptree.hpp>
#include <array>
#include <vector>
#include <utility>

namespace SpaceBallistics
{
  //=========================================================================//
  // "Ascent21b" Class:                                                      //
  //=========================================================================//
  class Ascent21b: public LVBase<Ascent21b>
  {
  public:
    //=======================================================================//
    // Consts:                                                               //
    //=======================================================================//
    // For Constrained Optimisation:
    // The max number of params used in Optimisation: 6:
    // [aAoAHat2, bAoAHat2, T3Rel, aAoAHat1, bAoAHat1, payLoadMassRel]:
    constexpr static int      NP             = 6;

    // Atmospheric Pressure at which Fairing Separation occurs (XXX: should it
    // be made configurable?):
    constexpr static Pressure FairingSepCond = Pressure(1.0);

    //=======================================================================//
    // Types:                                                                //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // The Base Class:                                                       //
    //-----------------------------------------------------------------------//
    using  Base = LVBase<Ascent21b>;
    friend class  LVBase<Ascent21b>;

    //-----------------------------------------------------------------------//
    // "FlightMode":                                                         //
    //-----------------------------------------------------------------------//
    // Soyuz-2.1b: parallel staging: the Stage1 Boosters and the Stage2 Core
    // burn together from the Lift-Off; Stage3 is ignited at the Stage2 Cut-
    // Off, prior to the Stage2 Separation. The modes are switched by the Stage
    // Models' times (in the Bwd integration, from Burn3 to Burn12):
    //
    enum class FlightMode: int
    {
      UNDEFINED = 0,
      Burn12    = 1,  // Stage1 + Stage2 (+ Stage3):     Lift-Off .. St1 Sep
      Burn2     = 2,  // Stage2          (+ Stage3):     St1 Sep  .. St2 CutOff
      Burn23    = 3,  // Stage2 (spent)  + Stage3 Burn:  St2 Cut  .. St2 Sep
      Burn3     = 4   // Stage3:                         St2 Sep  .. Insertion
    };

    static char const* ToString(FlightMode a_mode)
    {
      switch (a_mode)
      {
        case FlightMode::Burn12: return "Burn12";
        case FlightMode::Burn2 : return "Burn2";
        case FlightMode::Burn23: return "Burn23";
        case FlightMode::Burn3 : return "Burn3";
        default:                 return "UNDEFINED";
      }
    }

  private:
    //=======================================================================//
    // Data Flds:                                                            //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Over-All:                                                             //
    //-----------------------------------------------------------------------//
    // Non-Const Over-All Params (may be updated in the course of Optim-
    // isation):
    Mass                  m_payLoadMass;

    //-----------------------------------------------------------------------//
    // Mission Params:                                                       //
    //-----------------------------------------------------------------------//
    LenK                  m_Rins;        // Radius-Vector @ Orbital Insertion
    VelK                  m_Vins;        // LV Velocity   @ Orbital Insertion

    //-----------------------------------------------------------------------//
    // Flight Control Program Parameterisation:                              //
    //-----------------------------------------------------------------------//
    // AoA for Stage2: AoA(t) = t * (a * t - b),
    // where t <= 0 is our standard integration backward-running time    (t=0
    // corresponds  to the orbital insertion instant, and AoA=0 at that point);
    // NB: the (-b) coeff !!!
    Angle    const        m_maxAoA2;
    AngAcc                m_aAoA2;
    AngVel                m_bAoA2;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_aAoAHat2;
    double                m_bAoAHat2;

    // AoA for Stage1: for similarity with Stage2,
    // AoA(tau) = tau * (a * tau + b),
    // where tau  = t - tIgn1 >= 0, so AoA=0 @ tau=0 (Stage1 ignition):
    Angle    const        m_maxAoA1;
    AngAcc                m_aAoA1;
    AngVel                m_bAoA1;
    // Also, the Dim-Less params from which the above coeffs are computed:
    double                m_aAoAHat1;
    double                m_bAoAHat1;

    //-----------------------------------------------------------------------//
    // Transient Data (during flight path integration):                      //
    //-----------------------------------------------------------------------//
    FlightMode            m_mode;
    Time                  m_fairingSepTime; // Fairing Separation Time
    std::string           m_eventStr;
    Time                  m_nextOutputTime;

    // Vals which may be Constrained:
    Pressure              m_maxQ;           // Max Dynamic Pressure (Q)
    Pressure              m_sepQ;           // Q @ Stage1 separation
    double                m_maxLongG;       // Max Longitudinal G

    //-----------------------------------------------------------------------//
    // Soyuz-2.1b: the Stage3 BurnTime, the Segments, the Stage Models:      //
    //-----------------------------------------------------------------------//
    // The Closing Param, and the over-all flight time derived from it:
    Time                  m_T3;          // Stage3 BurnTime
    Time                  m_Tflight;     // Lift-Off .. Orbital Insertion

    // The AoA laws of "Ascent2": the Stage1 one applies to "Burn12" over the
    // domain 0 <= tau <= T1 ("Base::m_T1", set to Stage1SepTime); the Stage2
    // one, to the rest of the flight (Burn2, Burn23 and Burn3), over -T2 <= t
    // <= 0, T2 = Tflight - Stage1SepTime:
    Time                  m_T2;

    // The curr Segment (delimited by the Stage Models' events):
    Time                  m_segTop;         // Fwd Time of the curr Segment top
    Time                  m_segBot;         //   and of its bottom
    char const*           m_segEvent;       // Event @ "m_segTop" (or NULL)

    // The Stage Models' sum, cached for the last (Fwd Time, Pressure) pair: it
    // is used by "LVMass", "PropBurnRate" and "Thrust" at the same instant:
    struct Stack
    {
      Time                m_tau;
      Pressure            m_p;
      Density             m_rho;
      Mass                m_mass;       // Stages + InterStage + PayLoad
      MassRate            m_burnRate;   // >= 0
      ForceK              m_thrust;     // Along the LV axis
    };
    mutable Stack         m_stack;

  public:
    //=======================================================================//
    // Methods:                                                              //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    // Soyuz-2.1b: the Stage params come from the Stage Models, so only the Ctl
    // limits are left; the Launch Site is a "Soyuz2_LaunchPad":
    //
    Ascent21b
    (
      // Stage3:
      Time                    a_T3,              // Stage3 BurnTime
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

      // Integration / Output Params:
      Time                    a_ode_integr_step,
      std::ostream*           a_os,
      int                     a_log_level
    );

    //-----------------------------------------------------------------------//
    // Copy Ctor: Required for Optimisation:                                 //
    //-----------------------------------------------------------------------//
    Ascent21b(Ascent21b const&) = default;

    //-----------------------------------------------------------------------//
    // "SetCtlParams":                                                       //
    //-----------------------------------------------------------------------//
    // Sets the Ctl Params which are not set by the above Ctor (to avoid over-
    // complexity of the latter):
    //
    void SetCtlParams
    (
      double   a_aAoAHat2, double   a_bAoAHat2,
      double   a_aAoAHat1, double   a_bAoAHat1
    );

    //-----------------------------------------------------------------------//
    // "Run": Integrate the Ascent Trajectory:                               //
    //-----------------------------------------------------------------------//
    Base::RunRes Run();

  private:
    //=======================================================================//
    // Internal Helpers:                                                     //
    //=======================================================================//
    //-----------------------------------------------------------------------//
    // For ODE Integration:                                                  //
    //-----------------------------------------------------------------------//
    // "ODECB":
    // CallBack (invoked after the completion of each integration step):
    // Here FlightMode switching occurs, so this method is non-"const":
    //
    bool ODECB
    (
      StateV*       a_s,      Time a_t,      DStateV const& a_ds,
      StateV const& a_prev_s, Time a_prev_t, DStateV const& a_prev_ds
    );

    // "AeroDynForces":
    // Atmospheric Conditions and Aerodynamic Drag Force:
    //
    std::tuple<EAM::AtmConds, ForceK, ForceK>
    AeroDynForces(LenK a_r, VelK a_v, Angle a_AoA) const;

    // "PropBurnRate": (may be variable over time, >= 0):
    MassRate PropBurnRate(Mass a_m, LenK a_h, Time a_t) const;

    // "AoA": returns [Angle-of-Attack, ThrustVectorElevation]:
    std::pair<Angle,Angle> AoA(Time a_t, Angle a_psi) const;

    // "Thrust": Depends on the Mode, BurnRate and the Counter-Pressure:
    ForceK Thrust(MassRate a_burn_rate, Pressure a_p) const;

    // "LVMass": Current Mass (LV + PayLoad):
    Mass LVMass(StateV const& a_s, Time a_t) const;

    //-----------------------------------------------------------------------//
    // "FallingBackExn" Class:                                               //
    //-----------------------------------------------------------------------//
    // To be thrown if  the "ODECB"  has detected  that we are falling back to
    // Earth rather than ascending to orbit, so the integration process cannot
    // continue. It is an error cond, so it does not carry any info:
    //
    struct FallingBackExn
    {
    public:
      // Data Flds:
      LenK  const m_r;
      VelK  const m_Vr;
      VelK  const m_Vhor;
      Mass  const m_spentPropMass;
      Angle const m_phi;
      Time  const m_t;

      // Non-Default Ctor:
      FallingBackExn
      (
        LenK  a_r,   VelK a_Vr, VelK a_Vhor, Mass a_spent_prop_mass,
        Angle a_phi, Time a_t
      )
      : m_r             (a_r),
        m_Vr            (a_Vr),
        m_Vhor          (a_Vhor),
        m_spentPropMass (a_spent_prop_mass),
        m_phi           (a_phi),
        m_t             (a_t)
      {}
    };

    //-----------------------------------------------------------------------//
    // For Optimisation:                                                     //
    //-----------------------------------------------------------------------//
    // NB:   "ModifyLVParams" also calls "SetCtlParams", but is not limited  to
    // that: it also sets the Stage3 BurnTime and the PayLoadMass (the other
    // Mass and Thrust params are those of the Stage Models):
    //
    void ModifyLVParams
    (
      // Mass Params:
      Mass   a_payLoadMass,
      // Ctl Params:
      double a_aAoAHat2,    double a_bAoAHat2,
      Time   a_T3,
      double a_aAoAHat1,    double a_bAoAHat1
    );

    //-----------------------------------------------------------------------//
    // Soyuz-2.1b: the Vehicle from the Stage Models:                        //
    //-----------------------------------------------------------------------//
    // "SetAoACoeffs": The AoA coeffs from the Dim-Less params and "T{1,2}":
    void SetAoACoeffs();

    // "Tau": Fwd Flight Time for the Bwd Time "a_t", within the curr Segment:
    Time Tau(Time a_t) const;

    // "GetStack": The Stage Models' sum @ the Fwd Time "a_tau":
    Stack const& GetStack(Time a_tau, Pressure a_p, Density a_rho) const;

    //=======================================================================//
    // "NOMADEvaluator": Helper Class used in NOMAD Optimisation:            //
    //=======================================================================//
    class        NOMADEvaluator;
    friend class NOMADEvaluator;

  public:
    //=======================================================================//
    // "OptRes" Struct: The result of "FindOptimalAscentCtls" below:         //
    //=======================================================================//
    // It just contains the relevant flds taken from the main "Ascent21b" obj:
    //
    struct OptRes
    {
      //=====================================================================//
      // Data Flds:                                                          //
      //=====================================================================//
      //---------------------------------------------------------------------//
      // Stage2 Params:                                                      //
      //---------------------------------------------------------------------//
      // Stage2 AoA: Dim-Less Params:
      double    const m_aAoAHat2;         // Must be in [0 .. 1], default is 0
      double    const m_bAoAHat2;         // Must be in [0 .. 1], default is 0

      // Stage2 AoA: Dimensioned Coeffs: AoA2(t) = t * (a2 * t - b2), t <= 0,
      // t=0 is the orbital insertion time:
      AngAcc    const m_aAoA2;
      AngVel    const m_bAoA2;

      //---------------------------------------------------------------------//
      // Stage1 Params:                                                      //
      //---------------------------------------------------------------------//
      // Stage1 AoA: Dim-Less Params:
      double    const m_aAoAHat1;         // Must be in [0 .. 1], default is 0
      double    const m_bAoAHat1;         // Must be in [0 .. 1], default is 0

      // Stage1 AoA: Dimensioned Coeffs: AoA1(t) = tau * (a1 * tau + b1),
      // tau >= 0 is the time since Stage1 ignition:
      AngAcc    const m_aAoA1;
      AngVel    const m_bAoA1;

      //---------------------------------------------------------------------//
      // Over-All:                                                           //
      //---------------------------------------------------------------------//
      Mass      const m_payLoadMass;

      //---------------------------------------------------------------------//
      // Soyuz-2.1b: the AoA Domains and the Stage3 BurnTime:                //
      //---------------------------------------------------------------------//
      Time      const m_T2;               // AoA2: -T2 <= t   <= 0
      Time      const m_T1;               // AoA1:   0 <= tau <= T1
      Time      const m_T3;
      Time      const m_Tflight;
      Mass      const m_reserve3;         // Stage3 Propellant left

      // "FEASIBLE", "INFEASIBLE (NOT CLOSED)" (the least infeasible point,
      // see "Technical.WithInfeasibleResult") or "NOT OPTIMISED" (np=0):
      char const* const m_status;

      //=====================================================================//
      // Non-Default Ctor:                                                   //
      //=====================================================================//
      OptRes(Ascent21b const& a_asc, char const* a_status);

      //=====================================================================//
      // Output:                                                             //
      //=====================================================================//
      friend std::ostream& operator<<
            (std::ostream& a_os, OptRes const& a_res);
    };
    friend struct OptRes;

    //=======================================================================//
    // "FindOptimalAscentCtls":                                              //
    //=======================================================================//
    // Top-Level Optimisation Function.
    // Tries to find the Flight Control Params such that  the specified target
    // orbit is reached (with the orbital insertion occurring in the perigee,
    // if the orbit is elliptical) and that requires the min LV StartMass  (if
    // the PayloadMass is fixed)  or otherwise, achieves the max PayLoadMass.
    // Returns "OptRes" if successful, "nullopt" otherwise.
    // Can also perform the "final run" on the optimal params found, and return
    // the corresp "RunRes" (if specified in the Config.ini):
    //
    static std::pair<std::optional<OptRes>, std::optional<RunRes>>
    FindOptimalAscentCtls
    (
      // All are are given via the ConfigFile.ini, as there are quite a few of
      // them:
      std::string const& a_config_ini,
      std::ostream*      a_os   // May be NULL
    );

    // Same, with the Config already parsed (eg modified by the caller):
    static std::pair<std::optional<OptRes>, std::optional<RunRes>>
    FindOptimalAscentCtls
    (
      boost::property_tree::ptree const& a_config,
      std::ostream*                      a_os   // May be NULL
    );

    //=======================================================================//
    // Soyuz-2.1b: the Stage3 BurnTime Closure:                              //
    //=======================================================================//
    // The Stage Models run on the FORWARD flight time tau = t + Tflight, where
    // Tflight = Stage2CutOffTime + T3, and "T3" is the Stage3 (Block I) burn
    // time: at most the full burn of the Stage Model ("Stage3BurnDur"),  the
    // rest is the reserve. "T3" is the closing param of the trajectory:  the
    // Bwd integration must end at the Lift-Off time:
    //
    // The closure tolerance on the End Time Miss, used by "FindBurnTime3" and
    // the Optimiser:
    constexpr static Time     EndTimeTol     = 0.1_sec;

    //-----------------------------------------------------------------------//
    // "CloseRes": The result of "FindBurnTime3" below:                      //
    //-----------------------------------------------------------------------//
    struct CloseRes
    {
      bool      m_closed;   // Reached the Lift-Off State (see "FindBurnTime3")
      Time      m_T3;       // Stage3 Burn Time found
      Time      m_dT;       // End Time Miss: m_res.m_T + Tflight (0 is ideal)
      RunRes    m_res;      // Of the final run, at "m_T3"
      Mass      m_reserve3; // Stage3 Propellant left at its Cut-Off
      int       m_nRuns;    // Number of Runs performed
    };

    //-----------------------------------------------------------------------//
    // "SetBurnTime3", Accessors:                                            //
    //-----------------------------------------------------------------------//
    void SetBurnTime3(Time a_T3);

    Time GetBurnTime3    () const { return m_T3;             }
    Time GetFlightTime   () const { return m_Tflight;        }
    Time GetFairingSepTime() const { return m_fairingSepTime; }

    // Stage3 Propellant left (above the UnSpendable Remnants) at its Cut-Off:
    Mass GetReserve3() const;

    //-----------------------------------------------------------------------//
    // "FindBurnTime3":                                                      //
    //-----------------------------------------------------------------------//
    // 1-D search (bracketing and bisection) of the Stage3 BurnTime  such that
    // the Bwd integration ends at the Lift-Off time; the Ctl Params are kept.
    // "m_closed" is set if, in addition,  the end state is the Lift-Off one,
    // ie its Equivalent Velocity @ H=0 is within "a_max_start_v":
    //
    CloseRes FindBurnTime3(VelK a_max_start_v) const;
  };
}
// End namespace SpaceBallistics
