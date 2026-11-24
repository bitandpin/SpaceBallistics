// vim:ts=2:et
//===========================================================================//
//                    "Src/Missions/Ascent21b-Close.cpp":                    //
//    Ascent-to-Orbit of the "Soyuz-2.1b" LV: Closing by the Stage3 BurnTime //
//===========================================================================//
#include "SpaceBallistics/Missions/Ascent21b.h"
#include "SpaceBallistics/LVSC/Soyuz-2.1b/Stage3.h"
#include <algorithm>
#include <stdexcept>

namespace SpaceBallistics
{
//===========================================================================//
// "FindBurnTime3":                                                          //
//===========================================================================//
// The Bwd integration ends, depending on T3, either at the Singular Point
// (V=0) reached AFTER the Lift-Off time ("long" T3), or else ("short" T3)
// at the Lift-Off time with V > 0 (RunRC::FlameOut), or on the ground with
// V > 0 (RunRC::ZeroH). This classification is NOT monotonic in T3, so the
// "long" side is followed DOWN from the full burn in coarse steps, and the
// first "short" T3 met gives the bracket,  which is then bisected. The final
// run is on the "long" side of the bracket, ie it ends at the Lift-Off time
// up to the bisection tolerance; its "VT" is the remaining miss:
//
Ascent21b::CloseRes Ascent21b::FindBurnTime3(VelK a_max_start_v) const
{
  // XXX: Tolerances: on T3 (the bisection), and on the End Time Miss (the
  // log output interval, "EndTimeTol"); the number of coarse steps over the
  // T3 range:
  constexpr Time T3Tol  = 1e-3_sec;
  constexpr int  NSteps = 20;

  int nRuns = 0;

  // A single Run on a copy of this object (which must not have been Run):
  auto run =
    [this, &nRuns](Time a_T3, bool a_final) -> std::pair<RunRes, Time>
    {
      Ascent21b asc(*this);
      asc.m_logLevel = a_final ? m_logLevel : std::min(m_logLevel, 2);
      asc.SetBurnTime3(a_T3);
      ++nRuns;

      RunRes res;
      try
        { res = asc.Run(); }
      catch (std::exception const& exn)
      {
        if (m_os != nullptr && m_logLevel >= 1)
#         pragma omp critical(Output)
          *m_os << "# Ascent21b::FindBurnTime3: ERROR: T3="
                << a_T3.Magnitude() << " sec: " << exn.what() << std::endl;
      }
      if (m_os != nullptr && m_logLevel >= 2)
#       pragma omp critical(Output)
        *m_os << "# Ascent21b::FindBurnTime3: T3=" << a_T3.Magnitude()
              << " sec: RC="  << Base::ToString(res.m_rc) << ", T="
              << res.m_T.Magnitude()                 << " sec, Tflight="
              << asc.m_Tflight.Magnitude()           << " sec, VT="
              << res.m_VT.Magnitude()                << " km/sec"
              << std::endl;
      return std::make_pair(res, asc.m_Tflight);
    };

  // "Long" T3: the Singular Point is reached after the Lift-Off time (an Error
  // counts as "short"):
  auto isLong =
    [](std::pair<RunRes, Time> const& a_r) -> bool
    {
      return a_r.first.m_rc == RunRC::Singularity &&
             a_r.first.m_T  >  - a_r.second;
    };

  //-------------------------------------------------------------------------//
  // The Bracket:                                                            //
  //-------------------------------------------------------------------------//
  // NB: Within the shortest T3 allowed by "SetBurnTime3" (Stage3 burning alone
  // at the Insertion) and the full burn of the Stage Model:
  Time const tMin = Time(SC::Stage2SepTime) - Time(SC::Stage2CutOffTime) +
                    Base::m_odeIntegrStep;
  Time const step = (SC::Stage3BurnDur - tMin) / double(NSteps);
  Time       hi   = SC::Stage3BurnDur;
  Time       lo   = hi;

  if (isLong(run(hi, false)))
    while (true)
    {
      lo = std::max(hi - step, tMin);
      if (!isLong(run(lo, false)))
        break;
      hi = lo;
      if (hi <= tMin)
        break;       // Still "long" at the shortest T3
    }
  // Otherwise, even the full burn is "short", or an Error has occurred: the
  // final run is at "hi" anyway, to report the miss.

  //-------------------------------------------------------------------------//
  // Bisection:                                                              //
  //-------------------------------------------------------------------------//
  if (lo < hi)
    while (hi - lo > T3Tol)
    {
      Time mid  = 0.5 * (lo + hi);
      if (isLong(run(mid, false)))
        hi = mid;
      else
        lo = mid;
    }

  //-------------------------------------------------------------------------//
  // The Final Run (at the curr log level):                                  //
  //-------------------------------------------------------------------------//
  auto      fin = run(hi, true);
  Ascent21b asc(*this);
  asc.SetBurnTime3(hi);
  Time      dT  = fin.first.m_T + fin.second;

  bool closed =
    fin.first.m_rc != RunRC::Error && fin.first.m_VT <= a_max_start_v &&
    Abs(dT)        <= EndTimeTol;

  return CloseRes
  {
    .m_closed   = closed,
    .m_T3       = hi,
    .m_dT       = dT,
    .m_res      = fin.first,
    .m_reserve3 = asc.GetReserve3(),
    .m_nRuns    = nRuns
  };
}
}
// End namespace SpaceBallistics
