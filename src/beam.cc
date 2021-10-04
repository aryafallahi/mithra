/********************************************************************************************************
 *  undulator.cpp : Implementation of the functions for calculation of external field
 ********************************************************************************************************/

#include <string>
#include <fstream>

#include "solver.h"
#include "classes.h"

namespace MITHRA
{
  /* Calculate the fields of a static undulator.							*/
  void Solver::staticUndulator(UpdateBunchParallel& ubp, typename std::vector<Undulator>::iterator& iter)
  {
    if ( ubp.lz >= 0.0 && ubp.lz <= iter->length_ * iter->lu_ )
      {
	ubp.d1     = ub_.b0 * cosh(ub_.ku * ubp.ly) * sin(ub_.ku * ubp.lz) * gamma_;
	ubp.bt[0] += ubp.d1 * ub_.ct;
	ubp.bt[1] += ubp.d1 * ub_.st;
	ubp.bt[2] += ub_.b0 * sinh(ub_.ku * ubp.ly) * cos(ub_.ku * ubp.lz);

	ubp.d1    *= c0_ * beta_;
	ubp.et[1] +=   ubp.d1 * ub_.ct;
	ubp.et[0] += - ubp.d1 * ub_.st;
	ubp.et[2] += 0.0;
      }
    else if ( ubp.lz < 0.0 )
      {
	ubp.sz = exp( - pow( ub_.ku * ubp.lz , 2 ) / 2.0 );

	if ( iter != undulator_.begin() )
	  {
	    ubp.i  = iter - undulator_.begin() - 1;
	    ubp.r0 = undulator_[ubp.i].rb_ + undulator_[ubp.i].length_ * undulator_[ubp.i].lu_ - iter->rb_;
	    if ( ubp.lz < ubp.r0 || ubp.r0 == 0.0 ) ubp.sz = 0.0;
	    else
	      ubp.sz *= 0.35875 + 0.48829 * cos( PI * ubp.lz / ubp.r0 ) + 0.14128 * cos( 2.0 * PI * ubp.lz / ubp.r0 ) + 0.01168 * cos( 3.0 * PI * ubp.lz / ubp.r0 );
	  }

	ubp.d1     = ub_.b0 * cosh(ub_.ku * ubp.ly ) * ubp.sz * ub_.ku * ubp.lz * gamma_;
	ubp.bt[0] += ubp.d1 * ub_.ct;
	ubp.bt[1] += ubp.d1 * ub_.st;
	ubp.bt[2] += ub_.b0 * sinh(ub_.ku * ubp.ly) * ubp.sz;

	ubp.d1    *= c0_ * beta_;
	ubp.et[1] +=   ubp.d1 * ub_.ct;
	ubp.et[0] += - ubp.d1 * ub_.st;
	ubp.et[2] += 0.0;
      }
    else if ( ubp.lz > iter->length_ * iter->lu_ )
      {
	ubp.t0 = ubp.lz - iter->length_ * iter->lu_;

	ubp.sz = exp( - pow( ub_.ku *  ubp.t0 , 2 ) / 2.0 );

	if ( iter+1 != undulator_.end() )
	  {
	    ubp.i  = iter - undulator_.begin() + 1;
	    ubp.r0 = undulator_[ubp.i].rb_ - iter->rb_ - iter->length_ * iter->lu_;
	    if ( ubp.t0 > ubp.r0 || ubp.r0 == 0.0 ) ubp.sz = 0.0;
	    else
	      ubp.sz *= 0.35875 + 0.48829 * cos( PI * ubp.t0 / ubp.r0 ) + 0.14128 * cos( 2.0 * PI * ubp.t0 / ubp.r0 ) + 0.01168 * cos( 3.0 * PI * ubp.t0 / ubp.r0 );
	  }

	ubp.d1     = ub_.b0 * cosh(ub_.ku * ubp.ly ) * ubp.sz * ub_.ku * ubp.t0 * gamma_;
	ubp.bt[0] += ubp.d1 * ub_.ct;
	ubp.bt[1] += ubp.d1 * ub_.st;
	ubp.bt[2] += ub_.b0 * sinh(ub_.ku * ubp.ly ) * ubp.sz;

	ubp.d1    *= c0_ * beta_;
	ubp.et[1] +=   ubp.d1 * ub_.ct;
	ubp.et[0] += - ubp.d1 * ub_.st;
	ubp.et[2] += 0.0;
      }
  }

  /* Calculate the fields of a plane-wave.								*/
  template<class T>
  void Solver::planeWave(UpdateBunchParallel& ubp, T& s)
  {
    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal = s.signal_.self(ubp.tl, ubp.p0);

    /* Calculate the field only if the signal value is larger than a limit.                   		*/
    if ( fabs(ubp.tsignal) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Provide vector to store the electric field of the external field.                  		*/
    ubp.eT.mv( s.amplitude_ * ubp.tsignal, s.polarization_ );

    /* Provide vector to store the magnetic field of the external field.                  		*/
    ubp.bT = cross( s.direction_, s.polarization_ );
    ubp.bT.mv( s.amplitude_ * ubp.tsignal / c0_, ubp.bT );
  }

  /* Calculate the fields of a truncated plane-wave.							*/
  template<typename T>
  void Solver::planeWaveTruncated(UpdateBunchParallel& ubp, T& s)
  {
    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal = s.signal_.self(ubp.tl, ubp.p0);

    /* Calculate the field only if the signal value is larger than a limit.				*/
    if ( fabs(ubp.tsignal) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Calculate the transverse distance to the center line.   						*/
    ubp.x  = ubp.rv * s.polarization_;
    ubp.yv = cross( s.direction_, s.polarization_ );
    ubp.y  = ubp.rv * ubp.yv;

    /* Calculate the field only if the particle stays inside the truncated region.			*/
    if ( pow(ubp.x/s.radius_[0],2) + pow(ubp.y/s.radius_[1],2) > 1.0 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Provide vector to store the electric field of the undulator.                       		*/
    ubp.eT.mv( s.amplitude_ * ubp.tsignal, s.polarization_ );

    /* Provide vector to store the magnetic field of the undulator.                       		*/
    ubp.bT = cross( s.direction_, s.polarization_ );
    ubp.bT.mv( s.amplitude_ * ubp.tsignal / c0_, ubp.bT );
  }

  /* Calculate the fields of a gaussian beam.								*/
  template<typename T>
  void Solver::gaussianBeam(UpdateBunchParallel& ubp, T& s)
  {
    /* Calculate the transverse distance to the center line.                              		*/
    ubp.x  = ubp.rv * s.polarization_;
    ubp.yv = cross( s.direction_, s.polarization_ );
    ubp.y  = ubp.rv * ubp.yv;

    /* Calculate the relative radius of the beam.                 					*/
    ubp.wrp = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[0] * s.zR_[0] ) );
    ubp.wrs = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[1] * s.zR_[1] ) );

    ubp.x0 = ubp.x / ubp.wrp;
    ubp.y0 = ubp.y / ubp.wrs;

    /* Calculate the field only if the particle stays close enough to the center of the beam.		*/
    if ( fabs(ubp.x0) > 4.0 * s.radius_[0] ||  fabs(ubp.y0) > 4.0 * s.radius_[1] )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal = s.signal_.self(ubp.tl, ubp.p0);

    /* Calculate the field only if the signal value is larger than a limit.				*/
    if ( fabs(ubp.tsignal) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Get the required atan data.									*/
    ubp.atanP = atan( ubp.z / s.zR_[0] );
    ubp.atanS = atan( ubp.z / s.zR_[1] );

    /* Compute the transverse vector between the point and the reference point.           		*/
    ubp.p0  = 0.5 * ( ubp.atanP + ubp.atanS ) - PI * ubp.z / s.l_ * ( pow( ubp.x0 / s.zR_[0] , 2 ) + pow( ubp.y0 / s.zR_[1] , 2 ) );

    ubp.t   = exp( - pow( ubp.x0/s.radius_[0], 2) - pow( ubp.y0/s.radius_[1], 2) ) / sqrt(ubp.wrs*ubp.wrp);
    ubp.t  *= s.amplitude_;

    /* Retrieve signal value at corrected time.                                           		*/
    ubp.p1         = ubp.p0 - PI/2.0;
    ubp.tsignal    = s.signal_.self(ubp.tl, ubp.p1);
    ubp.ex.mv( ubp.t * ubp.tsignal,					s.polarization_ );
    ubp.by.mv( ubp.t * ubp.tsignal / c0_,				ubp.yv 		);

    ubp.p1         = ubp.p0 + ubp.atanP;
    ubp.tsignal    = s.signal_.self(ubp.tl, ubp.p1);
    ubp.ez.mv( ubp.t * ( - ubp.x0 / s.zR_[0] ) * ubp.tsignal, 		s.direction_);

    ubp.p1         = ubp.p0 + ubp.atanS;
    ubp.tsignal    = s.signal_.self(ubp.tl, ubp.p1);
    ubp.bz.mv( ubp.t * ( - ubp.y0 / s.zR_[1] ) / c0_ * ubp.tsignal, 	s.direction_);

    /* Calculate the total electric and magnetic field.                                   		*/
    ubp.eT = ubp.ex; ubp.eT += ubp.ez;
    ubp.bT = ubp.by; ubp.bT += ubp.bz;
  }

  /* Calculate the fields of a standing gaussian beam.							*/
  template<typename T>
  void Solver::superGaussianBeam(UpdateBunchParallel& ubp, T& s)
  {
    /* Calculate the transverse distance to the center line.                              		*/
    ubp.x  = ubp.rv * s.polarization_;
    ubp.yv = cross( s.direction_, s.polarization_ );
    ubp.y  = ubp.rv * ubp.yv;

    /* Calculate the relative radius of the beam.                 					*/
    ubp.wrp = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[0] * s.zR_[0] ) );
    ubp.wrs = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[1] * s.zR_[1] ) );

    /* Calculate the field only if the particle stays close enough to the center of the beam.		*/
    if ( ( fabs(ubp.x) - s.order_[0] * s.radius_[0] ) > 4.0 * s.radius_[0] * ubp.wrp ||
	 ( fabs(ubp.y) - s.order_[1] * s.radius_[1] ) > 4.0 * s.radius_[1] * ubp.wrs )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal  = s.signal_.self(ubp.tl,  ubp.p0);

    if ( fabs(ubp.tsignal) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Get the required atan data.									*/
    ubp.atanP = atan( ubp.z / s.zR_[0] );
    ubp.atanS = atan( ubp.z / s.zR_[1] );

    /*  Calculate the effective amplitude.								*/
    ubp.af = s.amplitude_ / sqrt(ubp.wrs*ubp.wrp);

    /* Initialize the values for E and B fields.							*/
    ubp.ex = 0.0; ubp.ez = 0.0; ubp.bz = 0.0; ubp.by = 0.0;

    /* Loop over elements of the super-gaussian beam and add their fields.				*/
    for ( int i = - s.order_[0]; i <= s.order_[0]; i++ )
      for ( int j = - s.order_[1]; j <= s.order_[1]; j++ )
	{
	  /* Compute the transverse vector between the point and the reference point.			*/
	  ubp.x0 = ( ubp.x - i * s.radius_[0] ) / ubp.wrp;
	  ubp.y0 = ( ubp.y - j * s.radius_[1] ) / ubp.wrs;

	  if ( fabs(ubp.x0) > 4.0 * s.radius_[0] ||  fabs(ubp.y0) > 4.0 * s.radius_[1] ) continue;

	  ubp.p0 = 0.5 * ( ubp.atanP + ubp.atanS ) - PI * ubp.z / s.l_ * ( pow( ubp.x0 / s.zR_[0] , 2 ) + pow( ubp.y / s.zR_[1] , 2 ) );
	  ubp.t  = ubp.af * exp( - pow( ubp.x0 / s.radius_[0], 2) - pow( ubp.y0 / s.radius_[1], 2) );

	  /* Retrieve signal value at corrected time.                                           	*/
	  ubp.p1         = ubp.p0 - PI/2.0;
	  ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);

	  ubp.ex.pmv( ubp.t * ubp.tsignal,					s.polarization_ );
	  ubp.by.pmv( ubp.t * ubp.tsignal / c0_,				ubp.yv 		);

	  ubp.p1         = ubp.p0 + ubp.atanP;
	  ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);

	  ubp.ez.pmv( ubp.t * ( - ubp.x0 / s.zR_[0] ) * ubp.tsignal,		s.direction_	);

	  ubp.p1         = ubp.p0 + ubp.atanS;
	  ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);

	  ubp.bz.pmv( ubp.t * ( - ubp.y0 / s.zR_[1] ) * ubp.tsignal / c0_,	s.direction_    );
	}

    /* Calculate the total electric and magnetic field.                                   		*/
    ubp.eT = ubp.ex; ubp.eT += ubp.ez;
    ubp.bT = ubp.by; ubp.bT += ubp.bz;
  }

  /* Calculate the fields of a standing plane wave.							*/
  template<typename T>
  void Solver::standingPlaneWave(UpdateBunchParallel& ubp, T& s)
  {
    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal  = s.signal_.self(ubp.tl,  ubp.p0);
    ubp.tsignalm = s.signal_.self(ubp.tlm, ubp.p0);
    ubp.tsignale = ubp.tsignal - ubp.tsignalm;
    ubp.tsignalb = ubp.tsignal + ubp.tsignalm;

    /* Calculate the fields only if the signal value is larger than a limit.				*/
    if ( fabs(ubp.tsignale) < 1.0e-6 && fabs(ubp.tsignalb) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    ubp.eT.mv( s.amplitude_ * ubp.tsignale, s.polarization_ );

    ubp.bT = cross( s.direction_, s.polarization_ );
    ubp.bT.mv( s.amplitude_ * ubp.tsignalb / c0_, ubp.bT );
  }

  /* Calculate the fields of a truncated standing plane wave.						*/
  template<typename T>
  void Solver::standingPlaneWaveTruncated(UpdateBunchParallel& ubp, T& s)
  {
    /* Calculate the transverse distance to the center line.   						*/
    ubp.x  = ubp.rv * s.polarization_;
    ubp.yv = cross( s.direction_, s.polarization_ );
    ubp.y  = ubp.rv * ubp.yv;

    if ( (ubp.x/s.radius_[0],2) + pow(ubp.y/s.radius_[1],2) > 1.0 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal  = s.signal_.self(ubp.tl,  ubp.p0);
    ubp.tsignalm = s.signal_.self(ubp.tlm, ubp.p0);
    ubp.tsignale = ubp.tsignal - ubp.tsignalm;
    ubp.tsignalb = ubp.tsignal + ubp.tsignalm;

    /* Calculate the fields only if the signal value is larger than a limit.				*/
    if ( fabs(ubp.tsignale) < 1.0e-6 && fabs(ubp.tsignalb) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    ubp.eT.mv( s.amplitude_ * ubp.tsignale, s.polarization_ );

    ubp.bT = cross( s.direction_, s.polarization_ );
    ubp.bT.mv( s.amplitude_ * ubp.tsignalb / c0_, ubp.bT );
  }

  /* Calculate the fields of a standing gaussian beam.							*/
  template<typename T>
  void Solver::standingGaussianBeam(UpdateBunchParallel& ubp, T& s)
  {
    /* Calculate the transverse distance to the center line.                              		*/
    ubp.x  = ubp.rv * s.polarization_;
    ubp.yv = cross( s.direction_, s.polarization_ );
    ubp.y  = ubp.rv * ubp.yv;

    /* Calculate the relative radius of the beam.                 					*/
    ubp.wrp = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[0] * s.zR_[0] ) );
    ubp.wrs = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[1] * s.zR_[1] ) );

    ubp.x0 = ubp.x / ubp.wrp;
    ubp.y0 = ubp.y / ubp.wrs;

    /* Calculate the field only if the particle stays close enough to the center of the beam.		*/
    if ( fabs(ubp.x0) > 4.0 * s.radius_[0] ||  fabs(ubp.y0) > 4.0 * s.radius_[1] )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal  = s.signal_.self(ubp.tl,  ubp.p0);
    ubp.tsignalm = s.signal_.self(ubp.tlm, ubp.p0);
    ubp.tsignale = ubp.tsignal - ubp.tsignalm;
    ubp.tsignalb = ubp.tsignal + ubp.tsignalm;

    if ( fabs(ubp.tsignale) < 1.0e-6 && fabs(ubp.tsignalb) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Get the required atan data.									*/
    ubp.atanP = atan( ubp.z / s.zR_[0] );
    ubp.atanS = atan( ubp.z / s.zR_[1] );

    /* Compute the transverse vector between the point and the reference point.           		*/
    ubp.p0  = 0.5 * ( ubp.atanP + ubp.atanS ) - PI * ubp.z / s.l_ * ( pow( ubp.x0 / s.zR_[0] , 2 ) + pow( ubp.y0 / s.zR_[1] , 2 ) );

    ubp.t   = exp( - pow( ubp.x0/s.radius_[0], 2) - pow( ubp.y0/s.radius_[1], 2) ) / sqrt(ubp.wrs*ubp.wrp);
    ubp.t  *= s.amplitude_;

    /* Retrieve signal value at corrected time.                                           		*/
    ubp.p1         = ubp.p0 - PI/2.0;
    ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);
    ubp.p1	   = ubp.p0 - PI/2.0;
    ubp.tsignalm   = s.signal_.self(ubp.tlm, ubp.p1);
    ubp.ex.mv( ubp.t * ( ubp.tsignal - ubp.tsignalm ),					s.polarization_ );
    ubp.by.mv( ubp.t / c0_ * ( ubp.tsignal + ubp.tsignalm ),				ubp.yv);

    ubp.p1         = ubp.p0 + ubp.atanP;
    ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);
    ubp.p1	   = -ubp.p1;
    ubp.tsignalm   = s.signal_.self(ubp.tlm, ubp.p1);
    ubp.ez.mv( ubp.t * ( - ubp.x0 / s.zR_[0] ) * ( ubp.tsignal - ubp.tsignalm ), 	s.direction_ );

    ubp.p1         = ubp.p0 + ubp.atanS;
    ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);
    ubp.p1	   = -ubp.p1;
    ubp.tsignalm   = s.signal_.self(ubp.tlm, ubp.p1);
    ubp.bz.mv( ubp.t * ( - ubp.y0 / s.zR_[1] ) / c0_ * ( ubp.tsignal - ubp.tsignalm ), 	s.direction_ );

    /* Calculate the total electric and magnetic field.                                   		*/
    ubp.eT = ubp.ex; ubp.eT += ubp.ez;
    ubp.bT = ubp.by; ubp.bT += ubp.bz;
  }

  /* Calculate the fields of a standing gaussian beam.							*/
  template<typename T>
  void Solver::standingSuperGaussianBeam(UpdateBunchParallel& ubp, T& s)
  {
    /* Calculate the transverse distance to the center line.                              		*/
    ubp.x  = ubp.rv * s.polarization_;
    ubp.yv = cross( s.direction_, s.polarization_ );
    ubp.y  = ubp.rv * ubp.yv;

    /* Calculate the relative radius of the beam.                 					*/
    ubp.wrp = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[0] * s.zR_[0] ) );
    ubp.wrs = sqrt( 1.0 + ubp.z * ubp.z / ( s.zR_[1] * s.zR_[1] ) );

    /* Calculate the field only if the particle stays close enough to the center of the beam.		*/
    if ( ( fabs(ubp.x) - s.order_[0] * s.radius_[0] ) > 4.0 * s.radius_[0] * ubp.wrp ||
	( fabs(ubp.y) - s.order_[1] * s.radius_[1] ) > 4.0 * s.radius_[1] * ubp.wrs )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Retrieve signal value at corrected time.                                               		*/
    ubp.tsignal  = s.signal_.self(ubp.tl,  ubp.p0);
    ubp.tsignalm = s.signal_.self(ubp.tlm, ubp.p0);
    ubp.tsignale = ubp.tsignal - ubp.tsignalm;
    ubp.tsignalb = ubp.tsignal + ubp.tsignalm;

    if ( fabs(ubp.tsignale) < 1.0e-6 && fabs(ubp.tsignalb) < 1.0e-6 )
      {
	ubp.eT = 0.0;
	ubp.bT = 0.0;
	return;
      }

    /* Get the required atan data.									*/
    ubp.atanP = atan( ubp.z / s.zR_[0] );
    ubp.atanS = atan( ubp.z / s.zR_[1] );

    /*  Calculate the effective amplitude.								*/
    ubp.af = s.amplitude_ / sqrt(ubp.wrs*ubp.wrp);

    /* Initialize the values for E and B fields.							*/
    ubp.ex = 0.0; ubp.ez = 0.0; ubp.bz = 0.0; ubp.by = 0.0;

    /* Loop over elements of the super-gaussian beam and add their fields.				*/
    for ( int i = - s.order_[0]; i <= s.order_[0]; i++ )
      for ( int j = - s.order_[1]; j <= s.order_[1]; j++ )
	{
	  /* Compute the transverse vector between the point and the reference point.			*/
	  ubp.x0 = ( ubp.x - i * s.radius_[0] ) / ubp.wrp;
	  ubp.y0 = ( ubp.y - j * s.radius_[1] ) / ubp.wrs;

	  if ( fabs(ubp.x0) > 4.0 * s.radius_[0] ||  fabs(ubp.y0) > 4.0 * s.radius_[1] ) continue;

	  ubp.p0 = 0.5 * ( ubp.atanP + ubp.atanS ) - PI * ubp.z / ubp.l * ( pow( ubp.x0 / s.zR_[0] , 2 ) + pow( ubp.y / s.zR_[1] , 2 ) );
	  ubp.t  = ubp.af * exp( - pow( ubp.x0 / s.radius_[0], 2) - pow( ubp.y0 / s.radius_[1], 2) );

	  /* Retrieve signal value at corrected time.                                           	*/
	  ubp.p1         = ubp.p0 - PI/2.0;
	  ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);
	  ubp.p1	 = ubp.p0 - PI/2.0;
	  ubp.tsignalm   = s.signal_.self(ubp.tlm, ubp.p1);

	  ubp.ex.pmv( ubp.t * ( ubp.tsignal - ubp.tsignalm ),					s.polarization_ 	);
	  ubp.by.pmv( ubp.t * ( ubp.tsignal + ubp.tsignalm ) / c0_,				ubp.yv 			);

	  ubp.p1         = ubp.p0 + ubp.atanP;
	  ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);
	  ubp.p1	 = -ubp.p1;
	  ubp.tsignalm   = s.signal_.self(ubp.tlm, ubp.p1);

	  ubp.ez.pmv( ubp.t * ( - ubp.x0 / s.zR_[0] ) * ( ubp.tsignal - ubp.tsignalm ),  	s.direction_    	);

	  ubp.p1         = ubp.p0 + ubp.atanS;
	  ubp.tsignal    = s.signal_.self(ubp.tl,  ubp.p1);
	  ubp.p1	 = -ubp.p1;
	  ubp.tsignalm   = s.signal_.self(ubp.tlm, ubp.p1);

	  ubp.bz.pmv( ubp.t * ( - ubp.y0 / s.zR_[1] ) * ( ubp.tsignal - ubp.tsignalm ) / c0_,	s.direction_    	);
	}

    /* Calculate the total electric and magnetic field.                                   		*/
    ubp.eT = ubp.ex; ubp.eT += ubp.ez;
    ubp.bT = ubp.by; ubp.bT += ubp.bz;
  }

}       /* End of namespace MITHRA.                                                                    	*/
