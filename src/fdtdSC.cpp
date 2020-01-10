/********************************************************************************************************
 *  fdtdSC.cpp : Implementation of the functions for real fdtd time marching solution class for the darius
 *  code with space-charge effect.
 ********************************************************************************************************/

#include <sys/time.h>

#include "fdtdSC.h"

namespace Darius
{
  FdTdSC::FdTdSC( Mesh& 				mesh,
		  Bunch& 				bunch,
		  Seed& 				seed,
		  std::vector<Undulator>&		undulator,
		  std::vector<ExtField>& 		extField,
		  std::vector<FreeElectronLaser>& 	FEL )
  : Solver ( mesh, bunch, seed, undulator, extField, FEL )
  {};

  /******************************************************************************************************
   * The function which is called for solving the fields in time domain.
   ******************************************************************************************************/

  void FdTdSC::solve ()
  {
    /* Declare the required variables in the calculations to avoid redundant data decalaration.		*/
    timeval           			simulationStart, simulationEnd;
    Double 				deltaTime, p = 0.0;
    std::stringstream 			printedMessage;
    std::vector<FieldVector<Double> >*	at;
    std::vector<Double>*		ft;
    std::list<Charge>::iterator 	iter;

    /* Transfer the whole quantities to the electron rest frame.					*/
    lorentzBoost();

    /* Before starting with the simulation the matrix for the field values and the coordinates should
     * be initialized. This is done based on the given mesh length and the mesh resolution in the mesh
     * structure.											*/
    initialize();

    /* Retrieve time when we start the electric update part                                  		*/
    gettimeofday(&simulationStart, NULL);

    /* Set the precision of number reports.								*/
    std::cout << std::fixed;
    std::cout << std::setprecision(3);

    /* Now update the field starting from time zero for the total simulation time.			*/
    printmessage(std::string(__FILE__), __LINE__, std::string("-> Run the time domain simulation ...") );
    while (time_ < mesh_.totalTime_)
      {
	/* Update the fields for one time step using the FDTD algorithm					*/
	fieldUpdate();

	/* If sampling of the field is enabled and the rhythm for sampling is achieved. Sample the
	 * field at the given position and save them into the file.					*/
	if ( seed_.sampling_ && fmod(time_, seed_.samplingRhythm_) < mesh_.timeStep_ && time_ > 0.0 ) fieldSample();

	/* If visualization of the field is enabled and the rhythm for visualization is achieved,
	 * visualize the fields and save the vtk data in the given file name.				*/
	for (unsigned int i = 0; i < seed_.vtk_.size(); i++)
	  {
	    if ( seed_.vtk_[i].sample_ && fmod(time_, seed_.vtk_[i].rhythm_) < mesh_.timeStep_ && time_ > 0.0 )
	      {
		if 		( seed_.vtk_[i].type_ == ALLDOMAIN ) fieldVisualizeAllDomain(i);
		else if 	( seed_.vtk_[i].type_ == INPLANE   ) fieldVisualizeInPlane(i);
	      }
	  }

	/* If profiling of the field is enabled and the time for profiling is achieved, write the field
	 * profile and save the data in the given file name.						*/
	if (seed_.profile_)
	  {
	    for (unsigned int i = 0; i < seed_.profileTime_.size(); i++)
	      if ( time_ - seed_.profileTime_[i] < mesh_.timeStep_ && time_ > seed_.profileTime_[i] )
		fieldProfile();

	    if ( fmod(time_, seed_.profileRhythm_) < mesh_.timeStep_ && time_ > 0.0 && seed_.profileRhythm_ != 0 )
	      fieldProfile();
	  }

	/* Reset the charge and current values to zero.							*/
	currentReset();

	/* For the sake of having correct charge conservation in the implementation of PIC model, we need
	 * to first update the charge motion with having the initial position saved in the memory. Then,
	 * the current and charge update should all happen using the very first and the very last charge
	 * positions. THIS IS VERY IMPORTANT AND SHOULD NOT BE CHANGED IN THE FUTURE.			*/

	/* Update the position and velocity parameters.							*/
	for (iter = iterQB_; iter != iterQE_; iter++)
	  {
	    iter->rnm  = iter->rnp;
	    iter->gbnm = iter->gbnp;
	  }

	/* Update the bunch till the time of the bunch properties reaches the time instant of the
	 * field.											*/
	for (Double t = 0.0; t < nUpdateBunch_; t += 1.0)
	  {
	    bunchUpdate();
	    timeBunch_ += bunch_.timeStep_;
	    ++nTimeBunch_;
	  }

	/* Update the values of the current.								*/
	currentUpdate();

	/* Communicate the current among processors.							*/
	currentCommunicate();

	/* If sampling of the bunch is enabled and the rhythm for sampling is achieved. Sample the
	 * bunch and save them into the file.								*/
	if ( bunch_.sampling_ && fmod(time_, bunch_.rhythm_) < mesh_.timeStep_ && time_ > 0.0 ) bunchSample();

	/* If visualization of the bunch is enabled and the rhythm for visualization is achieved,
	 * visualize the bunch and save the vtk data in the given file name.				*/
	if ( bunch_.bunchVTK_ && fmod(time_, bunch_.bunchVTKRhythm_) < mesh_.timeStep_ && time_ > 0.0 ) bunchVisualize();

	/* If profiling of the bunch is enabled and the time for profiling is achieved, write the bunch
	 * profile and save the data in the given file name.						*/
	if (bunch_.bunchProfile_ > 0)
	  {
	    for (unsigned int i = 0; i < (bunch_.bunchProfileTime_).size(); i++)
	      if ( time_ - bunch_.bunchProfileTime_[i] < mesh_.timeStep_ && time_ > bunch_.bunchProfileTime_[i] )
		bunchProfile();
	    if ( fmod(time_, bunch_.bunchProfileRhythm_) < mesh_.timeStep_ && time_ > 0.0 && bunch_.bunchProfileRhythm_ != 0.0 )
	      bunchProfile();
	  }

	/* If radiation power of the FEL output is enabled and the rhythm for sampling is achieved.
	 * Sample the radiation power at the given position and save them into the file.		*/
	powerSample(); powerVisualize();

	/* If radiation energy of the FEL output is enabled and the rhythm for sampling is achieved.
	 * Sample the radiation energy at the given position and save them into the file.		*/
	energySample();

	/* Shift the computed fields and the time points for the fields.				*/
	at    = anm1_;
	anm1_ = an_;
	an_   = anp1_;
	anp1_ = at;

	ft    = fnm1_;
	fnm1_ = fn_;
	fn_   = fnp1_;
	fnp1_ = ft;

	timem1_ += mesh_.timeStep_;
	time_   += mesh_.timeStep_;
	timep1_ += mesh_.timeStep_;
	++nTime_;

	gettimeofday(&simulationEnd, NULL);
	deltaTime  = ( simulationEnd.tv_usec - simulationStart.tv_usec ) / 1.0e6;
	deltaTime += ( simulationEnd.tv_sec - simulationStart.tv_sec );

	if ( rank_ == 0 && time_/mesh_.totalTime_ * 1000.0 > p )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Percentage of the simulation completed (%)      = ") +
			 stringify(time_/mesh_.totalTime_ * 100.0) );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Average calculation time for each time step (s) = ") +
			 stringify(deltaTime/(double)(nTime_))     );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Estimated remaining time (min)                  = ") +
			 stringify( (mesh_.totalTime_/time_ - 1) * deltaTime / 60 ) );
	    p += 1.0;
	  }
      }

    /* Finalize the calculations and the data saving.							*/
    finalize();
  }

  /******************************************************************************************************
   * Reset the currents to zero.
   ******************************************************************************************************/

  void FdTdSC::currentReset ()
  {
    Double*  jn = &jn_[0][0];
    Double*  je = &jn_[(long)N1N0_*np_-1][2];
    while ( jn != je )
      *(jn++) = 0.0;
    *je = 0.0;

    Double*  rn = &rn_[0];
    Double*  re = &rn_[(long)N1N0_*np_-1];
    while ( rn != re )
      *(rn++) = 0.0;
    *re = 0.0;
  }

  /******************************************************************************************************
   * Update the currents at cell points for the filed update.
   ******************************************************************************************************/

  void FdTdSC::currentUpdate ()
  {
    FieldVector<Double>*	        jn   = &jn_[0];
    Double*			        rn   = &rn_[0];
    bool                              bp, bm;
    std::list<Charge>::iterator       it = chargeVectorn_.begin();

    /* Now, a loop over the charges should be performed and the currents should be updated.           	*/
    for (it = chargeVectorn_.begin(); it != chargeVectorn_.end(); it++)
      {
	/* Find the current of the particle.                                                          	*/
	uc_.rp  = it->rnp;
	uc_.rm  = it->rnm;

	/* Save the flag detecting that the particle is in the domain of the processor.               	*/
	bp = ( uc_.rp[0] < xmax_ - uc_.dx && uc_.rp[0] > xmin_ + uc_.dx &&
	    uc_.rp[1] < ymax_ - uc_.dy && uc_.rp[1] > ymin_ + uc_.dy &&
	    uc_.rp[2] < zp_[1]         && uc_.rp[2] >= zp_[0] );
	bm = ( uc_.rm[0] < xmax_ - uc_.dx && uc_.rm[0] > xmin_ + uc_.dx &&
	    uc_.rm[1] < ymax_ - uc_.dy && uc_.rm[1] > ymin_ + uc_.dy &&
	    uc_.rm[2] < zp_[1]         && uc_.rm[2] >= zp_[0] );

	/* Continue the loop if none of the above conditions are met.                                 	*/
	if ( ! (bp || bm) ) continue;

	/* Get the charge of the particle.                                                            	*/
	uc_.q = it->q;

	/* Get the indices of the macro-particle in the next time step.                               	*/
	uc_.ip  = (int) floor( ( uc_.rp[0] - xmin_ ) / uc_.dx );
	uc_.jp  = (int) floor( ( uc_.rp[1] - ymin_ ) / uc_.dy );
	uc_.kp  = (int) floor( ( uc_.rp[2] - zmin_ ) / uc_.dz );

	/* Get the indices of the macro-particle in the previous time-step.                           	*/
	uc_.im  = (int) floor( ( uc_.rm[0] - xmin_ ) / uc_.dx );
	uc_.jm  = (int) floor( ( uc_.rm[1] - ymin_ ) / uc_.dy );
	uc_.km  = (int) floor( ( uc_.rm[2] - zmin_ ) / uc_.dz );

	/* Compute the relay point.                                                                   	*/
	uc_.r[0]  = std::min( std::min( uc_.im, uc_.ip ) * uc_.dx + uc_.dx + xmin_,
			      std::max( std::max( uc_.im, uc_.ip ) * uc_.dx + xmin_, 0.5 * (uc_.rm[0] + uc_.rp[0]) ) );
	uc_.r[1]  = std::min( std::min( uc_.jm, uc_.jp ) * uc_.dy + uc_.dy + ymin_,
			      std::max( std::max( uc_.jm, uc_.jp ) * uc_.dy + ymin_, 0.5 * (uc_.rm[1] + uc_.rp[1]) ) );
	uc_.r[2]  = std::min( std::min( uc_.km, uc_.kp ) * uc_.dz + uc_.dz + zmin_,
			      std::max( std::max( uc_.km, uc_.kp ) * uc_.dz + zmin_, 0.5 * (uc_.rm[2] + uc_.rp[2]) ) );

	/* Compute the charge fluxes.                                                                 	*/
	uc_.jcm  = uc_.r;
	uc_.jcm -= uc_.rm;
	uc_.jcp  = uc_.rp;
	uc_.jcp -= uc_.r;

	/* If the charge is outside the computational domain stop the simulation.                     	*/
	if ( bp )
	  {
	    /* Find the indices of the node whose current should be considered.                       	*/
	    uc_.m   = N1N0_ * ( uc_.kp - k0_ ) + N1_ * uc_.ip + uc_.jp;

	    uc_.dxp = modf( ( 0.5 * ( uc_.rp[0] + uc_.r[0] ) - xmin_ ) / uc_.dx , &uc_.c );
	    uc_.dyp = modf( ( 0.5 * ( uc_.rp[1] + uc_.r[1] ) - ymin_ ) / uc_.dy , &uc_.c );
	    uc_.dzp = modf( ( 0.5 * ( uc_.rp[2] + uc_.r[2] ) - zmin_ ) / uc_.dz , &uc_.c );

	    /* Calculate the contributions to the currents of each vertex.                            	*/
	    uc_.x1  = 1.0 - uc_.dxp;
	    uc_.x2  = uc_.dxp;
	    uc_.y1  = 1.0 - uc_.dyp;
	    uc_.y2  = uc_.dyp;
	    uc_.z1  = 1.0 - uc_.dzp;
	    uc_.z2  = uc_.dzp;

	    (*(jn+uc_.m)            )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+N1_)        )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+1  )        )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+N1_+1)      )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_)      )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_+N1_)  )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_+1)    )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_+N1_+1))[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcp[0];

	    (*(jn+uc_.m)            )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+N1_)        )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+1  )        )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+N1_+1)      )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_)      )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_+N1_)  )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_+1)    )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_+N1_+1))[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcp[1];

	    (*(jn+uc_.m)            )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+N1_)        )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+1  )        )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcp[2];
	    (*(jn+uc_.m+N1_+1)      )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_)      )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_+N1_)  )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_+1)    )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_+N1_+1))[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcp[2];

	    *(rn+uc_.m)                  += uc_.q * uc_.x1 * uc_.y1 * uc_.z1;
	    *(rn+uc_.m+N1_)        	   += uc_.q * uc_.x2 * uc_.y1 * uc_.z1;
	    *(rn+uc_.m+1  )              += uc_.q * uc_.x1 * uc_.y2 * uc_.z1;
	    *(rn+uc_.m+N1_+1)      	   += uc_.q * uc_.x2 * uc_.y2 * uc_.z1;
	    *(rn+uc_.m+N1N0_)            += uc_.q * uc_.x1 * uc_.y1 * uc_.z2;
	    *(rn+uc_.m+N1N0_+N1_)        += uc_.q * uc_.x2 * uc_.y1 * uc_.z2;
	    *(rn+uc_.m+N1N0_+1)          += uc_.q * uc_.x1 * uc_.y2 * uc_.z2;
	    *(rn+uc_.m+N1N0_+N1_+1)      += uc_.q * uc_.x2 * uc_.y2 * uc_.z2;
	  }

	/* If the charge is outside the computational domain stop the simulation.                     	*/
	if ( bm )
	  {
	    /* Find the indices of the node whose current should be considered.                       	*/
	    uc_.m   = N1N0_ * ( uc_.km - k0_ ) + N1_ * uc_.im + uc_.jm;

	    uc_.dxm = modf( ( 0.5 * ( uc_.rm[0] + uc_.r[0] ) - xmin_ ) / uc_.dx , &uc_.c );
	    uc_.dym = modf( ( 0.5 * ( uc_.rm[1] + uc_.r[1] ) - ymin_ ) / uc_.dy , &uc_.c );
	    uc_.dzm = modf( ( 0.5 * ( uc_.rm[2] + uc_.r[2] ) - zmin_ ) / uc_.dz , &uc_.c );

	    /* Calculate the contributions to the currents of each vertex.                            	*/
	    uc_.x1  = 1.0 - uc_.dxm;
	    uc_.x2  = uc_.dxm;
	    uc_.y1  = 1.0 - uc_.dym;
	    uc_.y2  = uc_.dym;
	    uc_.z1  = 1.0 - uc_.dzm;
	    uc_.z2  = uc_.dzm;

	    (*(jn+uc_.m)            )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+N1_)        )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+1  )        )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+N1_+1)      )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_)      )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_+N1_)  )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_+1)    )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_+N1_+1))[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcm[0];

	    (*(jn+uc_.m)            )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+N1_)        )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+1  )        )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+N1_+1)      )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_)      )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_+N1_)  )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_+1)    )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_+N1_+1))[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcm[1];

	    (*(jn+uc_.m)            )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+N1_)        )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+1  )        )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcm[2];
	    (*(jn+uc_.m+N1_+1)      )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_)      )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_+N1_)  )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_+1)    )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_+N1_+1))[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcm[2];

	    *(rn+uc_.m)                  += uc_.q * uc_.x1 * uc_.y1 * uc_.z1;
	    *(rn+uc_.m+N1_)              += uc_.q * uc_.x2 * uc_.y1 * uc_.z1;
	    *(rn+uc_.m+1  )              += uc_.q * uc_.x1 * uc_.y2 * uc_.z1;
	    *(rn+uc_.m+N1_+1)            += uc_.q * uc_.x2 * uc_.y2 * uc_.z1;
	    *(rn+uc_.m+N1N0_)            += uc_.q * uc_.x1 * uc_.y1 * uc_.z2;
	    *(rn+uc_.m+N1N0_+N1_)        += uc_.q * uc_.x2 * uc_.y1 * uc_.z2;
	    *(rn+uc_.m+N1N0_+1)          += uc_.q * uc_.x1 * uc_.y2 * uc_.z2;
	    *(rn+uc_.m+N1N0_+N1_+1)      += uc_.q * uc_.x2 * uc_.y2 * uc_.z2;
	  }
      }
  }

  /******************************************************************************************************
   * Communicate the currents among different processors.
   ******************************************************************************************************/

  void FdTdSC::currentCommunicate ()
  {
    int                               msgtag9 = 9, msgtag10 = 10;
    MPI_Status                        status;
    std::list<Charge>::iterator       it = chargeVectorn_.begin();

    /* Add the contribution of each processor to the charge and current density at the boundaries.    	*/
    if (rank_ != 0)
      {
	MPI_Send(&jn_[0][0],          3*N1N0_,MPI_DOUBLE,rank_-1,msgtag9, MPI_COMM_WORLD);
	MPI_Send(&rn_[0],               N1N0_,MPI_DOUBLE,rank_-1,msgtag10,MPI_COMM_WORLD);
      }

    if (rank_ != size_ - 1)
      {
	MPI_Recv(&uc_.jt[0][0],       3*N1N0_,MPI_DOUBLE,rank_+1,msgtag9, MPI_COMM_WORLD,&status);
	MPI_Recv(&uc_.rt[0],            N1N0_,MPI_DOUBLE,rank_+1,msgtag10,MPI_COMM_WORLD,&status);

	for (int i = 0; i < N1N0_; i++)
	  {
	    jn_[(np_-2)*N1N0_+i] += uc_.jt[i];
	    rn_[(np_-2)*N1N0_+i] += uc_.rt[i];
	  }
      }

    /* Now that the charge and current densities are deposited, remove the out of domain charges from the
     * list.                                                                                      	*/
    it = chargeVectorn_.begin();
    while ( it != chargeVectorn_.end() )
      {
	if      ( it->rnp[2] <  zp_[0] && rank_ != 0 )
	  it = chargeVectorn_.erase(it);
	else if ( it->rnp[2] >= zp_[1] && rank_ != size_ - 1 )
	  it = chargeVectorn_.erase(it);
	else
	  ++it;
      }

    /* Initialize the end and begin of the charge vector iterator.                                    	*/
    iterQB_ = chargeVectorn_.begin();
    iterQE_ = chargeVectorn_.end();
  }

  /******************************************************************************************************
   * Update the fields for one time-step
   ******************************************************************************************************/

  void FdTdSC::fieldUpdate ()
  {
    /* Define the values temporally needed for updating the fields.					*/
    unsigned int 		i, j, k;
    long int                  m, l;
    MPI_Status 		status;
    int 			msgtag1 = 1, msgtag2 = 2, msgtag3 = 3, msgtag4 = 4;
    int                       msgtag5 = 5, msgtag6 = 6, msgtag7 = 7, msgtag8 = 8;
    FieldVector<Double>	atemp; atemp = 0.0;

    uf_.anp1 = &(*anp1_)[0][0];
    uf_.an   = &(*an_)  [0][0];
    uf_.anm1 = &(*anm1_)[0][0];
    uf_.jn   = &jn_[0][0];
    uf_.en   = &en_[0][0];
    uf_.bn   = &bn_[0][0];

    uf_.fnp1 = &(*fnp1_)[0];
    uf_.fn   = &(*fn_)  [0];
    uf_.fnm1 = &(*fnm1_)[0];
    uf_.rn   = &rn_[0];

    const long int M0  = N1_;
    const long int M1  = N1N0_;
    const long int M2  = N1_+N1N0_;
    const long int M3  = N1_-N1N0_;
    const long int M4  = 1+N1N0_;
    const long int M5  = 1-N1N0_;
    const long int M6  = 1+N1_;
    const long int M7  = 1-N1_;
    const long int M8  = 1+N1_+N1N0_;
    const long int M9  = 1+N1_-N1N0_;
    const long int M10 = 1-N1_+N1N0_;
    const long int M11 = 1-N1_-N1N0_;

    const long int L0  = 3*N1_;
    const long int L1  = 3*N1N0_;
    const long int L2  = 3*(N1_+N1N0_);
    const long int L3  = 3*(N1_-N1N0_);
    const long int L4  = 3*(1+N1N0_);
    const long int L5  = 3*(1-N1N0_);
    const long int L6  = 3*(1+N1_);
    const long int L7  = 3*(1-N1_);
    const long int L8  = 3*(1+N1_+N1N0_);
    const long int L9  = 3*(1+N1_-N1N0_);
    const long int L10 = 3*(1-N1_+N1N0_);
    const long int L11 = 3*(1-N1_-N1N0_);

    /* First set all the particle-in-cell flags equal to zero.						*/
    std::vector<bool>::iterator  pn = pic_.begin();
    std::vector<bool>::iterator  pe = pic_.end();
    while ( pn != pe ) *(pn++) = false;

    /* Loop over the points in the mesh and update the fields by one time step. Since, the time update
     * for the points in the computational domain are different from the ones on the boundary, we will
     * do a loop first on the internal points. After all of them are updated, the points on the
     * boundary will be updated accordingly.								*/
    if ( mesh_.solver_ == NSFD )
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    for (unsigned k = 1; k < uf_.npm1; k++)
	      {
		m = N1N0_ * k + N1_ * i + j;
		l = 3 * m;

		uf_.af.advanceMagneticPotentialNSFD(
		    uf_.anp1+l,    uf_.anm1+l,    uf_.an+l,
		    uf_.an  +l+L0, uf_.an  +l+L2, uf_.an+l+L3,
		    uf_.an  +l-L0, uf_.an  +l-L3, uf_.an+l-L2,
		    uf_.an  +l+3 , uf_.an  +l+L4, uf_.an+l+L5,
		    uf_.an  +l-3 , uf_.an  +l-L5, uf_.an+l-L4,
		    uf_.an  +l+L1, uf_.an  +l-L1, uf_.jn+l);

		uf_.af.advanceScalarPotentialNSFD(
		    uf_.fnp1+m,    uf_.fnm1+m,    uf_.fn+m,
		    uf_.fn  +m+M0, uf_.fn  +m+M2, uf_.fn+m+M3,
		    uf_.fn  +m-M0, uf_.fn  +m-M3, uf_.fn+m-M2,
		    uf_.fn  +m+1 , uf_.fn  +m+M4, uf_.fn+m+M5,
		    uf_.fn  +m-1 , uf_.fn  +m-M5, uf_.fn+m-M4,
		    uf_.fn  +m+M1, uf_.fn  +m-M1, uf_.rn+m);
	      }
      }
    else if ( mesh_.solver_ == FD )
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    for (unsigned k = 1; k < uf_.npm1; k++)
	      {
		m = N1N0_ * k + N1_ * i + j;
		l = 3 * m;

		uf_.af.advanceMagneticPotentialFD(
		    uf_.anp1+l,    uf_.anm1+l,    uf_.an+l,
		    uf_.an  +l+L0, uf_.an  +l+L2, uf_.an+l+L3,
		    uf_.an  +l-L0, uf_.an  +l-L3, uf_.an+l-L2,
		    uf_.an  +l+3 , uf_.an  +l+L4, uf_.an+l+L5,
		    uf_.an  +l-3 , uf_.an  +l-L5, uf_.an+l-L4,
		    uf_.an  +l+L1, uf_.an  +l-L1, uf_.jn+l);

		uf_.af.advanceScalarPotentialFD(
		    uf_.fnp1+m,    uf_.fnm1+m,    uf_.fn+m,
		    uf_.fn  +m+M0, uf_.fn  +m+M2, uf_.fn+m+M3,
		    uf_.fn  +m-M0, uf_.fn  +m-M3, uf_.fn+m-M2,
		    uf_.fn  +m+1 , uf_.fn  +m+M4, uf_.fn+m+M5,
		    uf_.fn  +m-1 , uf_.fn  +m-M5, uf_.fn+m-M4,
		    uf_.fn  +m+M1, uf_.fn  +m-M1, uf_.rn+m);
	      }
      }

    /* If the amplitude of the seed exceeds a certain limit inject the seed into the computational domain
     * through TF/SF boundaries.									*/
    if ( fabs(seed_.amplitude_) > 1.0e-50 )
      {
	unsigned int KI = ( rank_ == 0         ) ? 2       : 1;
	unsigned int KF = ( rank_ == size_ - 1 ) ? np_ - 2 : np_ - 1;

	for (int j = 2; j < N1_-2; j++)
	  for (unsigned k = KI; k < KF; k++)
	    {
	      i = 1;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m+N1_],	time_, atemp); (*anp1_)[m].mmv(uf_.a[1], atemp);
	      i = 2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m-N1_],	time_, atemp); (*anp1_)[m].pmv(uf_.a[1], atemp);
	      i = N0_-2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m-N1_],	time_, atemp); (*anp1_)[m].mmv(uf_.a[1], atemp);
	      i = N0_-3;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m+N1_],	time_, atemp); (*anp1_)[m].pmv(uf_.a[1], atemp);
	    }

	for (int i = 2; i < N0_-2; i++)
	  for (unsigned k = KI; k < KF; k++)
	    {
	      j = 1;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m+1], time_, atemp); (*anp1_)[m].mmv(uf_.a[2], atemp);
	      j = 2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m-1], time_, atemp); (*anp1_)[m].pmv(uf_.a[2], atemp);
	      j = N1_-2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m-1], time_, atemp); (*anp1_)[m].mmv(uf_.a[2], atemp);
	      j = N1_-3;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(r_[m+1],		time_, atemp); (*anp1_)[m].pmv(uf_.a[2], atemp);
	    }

	if ( rank_ == 0 )
	  {
	    for (int i = 2; i < N0_-2; i++)
	      for (int j = 2; j < N1_-2; j++)
		{
		  k = 1;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(r_[m+N1N0_], time_, atemp); (*anp1_)[m].mmv(uf_.a[3], atemp);
		  k = 2;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(r_[m-N1N0_], time_, atemp); (*anp1_)[m].pmv(uf_.a[3], atemp);
		}
	  }

	if ( rank_ == size_-1 )
	  {
	    for (int i = 2; i < N0_-2; i++)
	      for (int j = 2; j < N1_-2; j++)
		{
		  k = np_-2;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(r_[m-N1N0_],	time_, atemp); (*anp1_)[m].mmv(uf_.a[3], atemp);
		  k = np_-3;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(r_[m+N1N0_],	time_, atemp); (*anp1_)[m].pmv(uf_.a[3], atemp);
		}
	  }
      }

    /* Loop over the points in the mesh on the x = xmin boundary and update the fields using the first
     * order absorbing boundary condition.								*/
    uf_.af.ufB_ = &uf_.bB[0];
    for (unsigned j = 1; j < uf_.N1m1; j++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  m = N1N0_ * k + j;
	  l = 3 * m;

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l+L0,	uf_.an  +l+L0, 	uf_.anp1+l+L0,
	      uf_.an  +l+L6,  uf_.an  +l-L7,	uf_.an  +l+L2,
	      uf_.an  +l+L3, 	uf_.an  +l+3,   uf_.an  +l-3,
	      uf_.an  +l+L1,  uf_.an  +l-L1);
	  uf_.af.advanceBoundaryS(
	      uf_.fnp1+m,	uf_.fnm1+m,	uf_.fn  +m,
	      uf_.fnm1+m+M0,	uf_.fn  +m+M0, 	uf_.fnp1+m+M0,
	      uf_.fn  +m+M6,  uf_.fn  +m-M7,	uf_.fn  +m+M2,
	      uf_.fn  +m+M3, 	uf_.fn  +m+1,   uf_.fn  +m-1,
	      uf_.fn  +m+M1,  uf_.fn  +m-M1);
	}

    /* Loop over the points in the mesh on the x = xmax boundary and update the fields using the first
     * order absorbing boundary condition.								*/
    for (unsigned j = 1; j < uf_.N1m1; j++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  m = N1N0_ * k + N1N0_ - N1_ + j;
	  l = 3 * m;

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l-L0,	uf_.an  +l-L0,	uf_.anp1+l-L0,
	      uf_.an  +l+L7,	uf_.an  +l-L6,	uf_.an  +l-L3,
	      uf_.an  +l-L2, 	uf_.an  +l+3,	uf_.an  +l-3,
	      uf_.an  +l+L1,  uf_.an  +l-L1 );
	  uf_.af.advanceBoundaryS(
	      uf_.fnp1+m,	uf_.fnm1+m,	uf_.fn  +m,
	      uf_.fnm1+m-M0,	uf_.fn  +m-M0,	uf_.fnp1+m-M0,
	      uf_.fn  +m+M7,	uf_.fn  +m-M6,	uf_.fn  +m-M3,
	      uf_.fn  +m-M2, 	uf_.fn  +m+1,	uf_.fn  +m-1,
	      uf_.fn  +m+M1,  uf_.fn  +m-M1 );
	}

    /* Loop over the points in the mesh on the y = ymin boundary and update the fields using the first
     * order absorbing boundary condition.								*/
    uf_.af.ufB_ = &uf_.cB[0];
    for (unsigned i = 1; i < uf_.N0m1; i++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  m = N1N0_ * k + N1_* i;
	  l = 3 * m;

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l+3,   uf_.an  +l+3,   uf_.anp1+l+3,
	      uf_.an  +l+L6,  uf_.an  +l+L7, 	uf_.an  +l+L4,
	      uf_.an  +l+L5, 	uf_.an  +l+L0,  uf_.an  +l-L0,
	      uf_.an  +l+L1,  uf_.an  +l-L1);
	  uf_.af.advanceBoundaryS(
	      uf_.fnp1+m,	uf_.fnm1+m,	uf_.fn  +m,
	      uf_.fnm1+m+1,   uf_.fn  +m+1,   uf_.fnp1+m+1,
	      uf_.fn  +m+M6,  uf_.fn  +m+M7, 	uf_.fn  +m+M4,
	      uf_.fn  +m+M5, 	uf_.fn  +m+M0,  uf_.fn  +m-M0,
	      uf_.fn  +m+M1,  uf_.fn  +m-M1);
	}

    /* Loop over the points in the mesh on the y = ymax boundary and update the fields using the first
     * order absorbing boundary condition.								*/
    for (unsigned i = 1; i < uf_.N0m1; i++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  m = N1N0_ * k + N1_* i + N1_ - 1;
	  l = 3 * m;

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l-3, 	uf_.an  +l-3,  	uf_.anp1+l-3,
	      uf_.an  +l-L7,  uf_.an  +l-L6, 	uf_.an  +l-L5,
	      uf_.an  +l-L4,	uf_.an  +l+L0,  uf_.an  +l-L0,
	      uf_.an  +l+L1,  uf_.an  +l-L1);
	  uf_.af.advanceBoundaryS(
	      uf_.fnp1+m,	uf_.fnm1+m,	uf_.fn  +m,
	      uf_.fnm1+m-1, 	uf_.fn  +m-1,  	uf_.fnp1+m-1,
	      uf_.fn  +m-M7,  uf_.fn  +m-M6, 	uf_.fn  +m-M5,
	      uf_.fn  +m-M4,	uf_.fn  +m+M0,  uf_.fn  +m-M0,
	      uf_.fn  +m+M1,  uf_.fn  +m-M1);
	}

    /* Loop over the points in the mesh on the z = zmin boundary and update the fields using the first
     * order absorbing boundary condition.								*/
    uf_.af.ufB_ = &uf_.dB[0];
    if ( rank_ == 0 )
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    {
	      m = N1_ * i + j;
	      l = 3 * m;

	      uf_.af.advanceBoundaryF(
		  uf_.anp1+l,		uf_.anm1+l,		uf_.an  +l,
		  uf_.anm1+l+L1,     	uf_.an  +l+L1,   	uf_.anp1+l+L1,
		  uf_.an  +l+L2,   	uf_.an  +l-L3, 		uf_.an  +l+L4,
		  uf_.an  +l-L5, 	uf_.an  +l+L0,     	uf_.an  +l-L0,
		  uf_.an  +l+3,     	uf_.an  +l-3);
	      uf_.af.advanceBoundaryS(
		  uf_.fnp1+m,		uf_.fnm1+m,		uf_.fn  +m,
		  uf_.fnm1+m+M1,     	uf_.fn  +m+M1,   	uf_.fnp1+m+M1,
		  uf_.fn  +m+M2,   	uf_.fn  +m-M3, 		uf_.fn  +m+M4,
		  uf_.fn  +m-M5, 	uf_.fn  +m+M0,     	uf_.fn  +m-M0,
		  uf_.fn  +m+1,     	uf_.fn  +m-1);
	    }
      }

    /* Loop over the points in the mesh on the z = zmax boundary and update the fields using the first
     * order absorbing boundary condition.								*/
    if (rank_ == size_ - 1)
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    {
	      m = N1N0_ * uf_.npm1 + N1_ * i + j;
	      l = 3 * m;

	      uf_.af.advanceBoundaryF(
		  uf_.anp1+l,	uf_.anm1+l,		uf_.an  +l,
		  uf_.anm1+l-L1,     uf_.an  +l-L1,   	uf_.anp1+l-L1,
		  uf_.an  +l+L3,   	uf_.an  +l-L2, 		uf_.an  +l+L5,
		  uf_.an  +l-L4, 	uf_.an  +l+L0,     	uf_.an  +l-L0,
		  uf_.an  +l+3,    	uf_.an  +l-3);
	      uf_.af.advanceBoundaryS(
		  uf_.fnp1+m,	uf_.fnm1+m,		uf_.fn  +m,
		  uf_.fnm1+m-M1,     uf_.fn  +m-M1,   	uf_.fnp1+m-M1,
		  uf_.fn  +m+M3,   	uf_.fn  +m-M2, 		uf_.fn  +m+M5,
		  uf_.fn  +m-M4, 	uf_.fn  +m+M0,     	uf_.fn  +m-M0,
		  uf_.fn  +m+1,    	uf_.fn  +m-1);
	    }
      }

    if ( mesh_.truncationOrder_ == 2 )
      {

	/* Loop over the edge points in the mesh on the x = (xmin,xmax) and y = (ymin,ymax) boundary and
	 * update the fields using the first order absorbing boundary condition. To understand the
	 * following lines of codes, it is better to list the values of Li at the side.			*/
	uf_.af.ufB_ = &uf_.eE[0];
	for (unsigned k = 1; k < uf_.npm1; k++)
	  {
	    m = N1N0_ * k;
	    l = 3 * m;

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l+L0,	uf_.an+l+L0,	uf_.anm1+l+L0,
		uf_.anp1+l+3,		uf_.an+l+3,	uf_.anm1+l+3,
		uf_.anp1+l+L6,	uf_.an+l+L6,	uf_.anm1+l+L6,
		uf_.an  +l-L1,	uf_.an+l+L3,	uf_.an  +l+L5,	uf_.an+l+L9,
		uf_.an  +l+L1,	uf_.an+l+L2,	uf_.an  +l+L4,	uf_.an+l+L8);
	    uf_.af.advanceEdgeS(
		uf_.fnp1+m,		uf_.fn+m,	uf_.fnm1+m,
		uf_.fnp1+m+M0,	uf_.fn+m+M0,	uf_.fnm1+m+M0,
		uf_.fnp1+m+1,		uf_.fn+m+1,	uf_.fnm1+m+1,
		uf_.fnp1+m+M6,	uf_.fn+m+M6,	uf_.fnm1+m+M6,
		uf_.fn  +m-M1,	uf_.fn+m+M3,	uf_.fn  +m+M5,	uf_.fn+m+M9,
		uf_.fn  +m+M1,	uf_.fn+m+M2,	uf_.fn  +m+M4,	uf_.fn+m+M8);

	    m = N1N0_ * k + N1_ * uf_.N0m1;
	    l = 3 * m;

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l-L0,	uf_.an+l-L0,	uf_.anm1+l-L0,
		uf_.anp1+l+3,		uf_.an+l+3,	uf_.anm1+l+3,
		uf_.anp1+l+L7,	uf_.an+l+L7,	uf_.anm1+l+L7,
		uf_.an  +l-L1,	uf_.an+l-L2,	uf_.an  +l+L5,	uf_.an+l+L11,
		uf_.an  +l+L1,	uf_.an+l-L3,	uf_.an  +l+L4,	uf_.an+l+L10);
	    uf_.af.advanceEdgeS(
		uf_.fnp1+m,		uf_.fn+m,	uf_.fnm1+m,
		uf_.fnp1+m-M0,	uf_.fn+m-M0,	uf_.fnm1+m-M0,
		uf_.fnp1+m+1,		uf_.fn+m+1,	uf_.fnm1+m+1,
		uf_.fnp1+m+M7,	uf_.fn+m+M7,	uf_.fnm1+m+M7,
		uf_.fn  +m-M1,	uf_.fn+m-M2,	uf_.fn  +m+M5,	uf_.fn+m+M11,
		uf_.fn  +m+M1,	uf_.fn+m-M3,	uf_.fn  +m+M4,	uf_.fn+m+M10);

	    m = N1N0_ * k + uf_.N1m1;
	    l = 3 * m;

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l+L0,	uf_.an+l+L0,	uf_.anm1+l+L0,
		uf_.anp1+l-3,		uf_.an+l-3,	uf_.anm1+l-3,
		uf_.anp1+l-L7,	uf_.an+l-L7,	uf_.anm1+l-L7,
		uf_.an  +l-L1,	uf_.an+l+L3,	uf_.an  +l-L4,	uf_.an+l-L10,
		uf_.an  +l+L1,	uf_.an+l+L2,	uf_.an  +l-L5,	uf_.an+l-L11);
	    uf_.af.advanceEdgeS(
		uf_.fnp1+m,		uf_.fn+m,	uf_.fnm1+m,
		uf_.fnp1+m+M0,	uf_.fn+m+M0,	uf_.fnm1+m+M0,
		uf_.fnp1+m-1,		uf_.fn+m-1,	uf_.fnm1+m-1,
		uf_.fnp1+m-M7,	uf_.fn+m-M7,	uf_.fnm1+m-M7,
		uf_.fn  +m-M1,	uf_.fn+m+M3,	uf_.fn  +m-M4,	uf_.fn+m-M10,
		uf_.fn  +m+M1,	uf_.fn+m+M2,	uf_.fn  +m-M5,	uf_.fn+m-M11);

	    m = N1N0_ * k + N1_ * uf_.N0m1 + uf_.N1m1;
	    l = 3 * m;

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l-L0,	uf_.an+l-L0,	uf_.anm1+l-L0,
		uf_.anp1+l-3,		uf_.an+l-3,	uf_.anm1+l-3,
		uf_.anp1+l-L6,	uf_.an+l-L6,	uf_.anm1+l-L6,
		uf_.an  +l-L1,	uf_.an+l-L2,	uf_.an  +l-L4,	uf_.an+l-L8,
		uf_.an  +l+L1,	uf_.an+l-L3,	uf_.an  +l-L5,	uf_.an+l-L9);
	    uf_.af.advanceEdgeS(
		uf_.fnp1+m,		uf_.fn+m,	uf_.fnm1+m,
		uf_.fnp1+m-M0,	uf_.fn+m-M0,	uf_.fnm1+m-M0,
		uf_.fnp1+m-1,		uf_.fn+m-1,	uf_.fnm1+m-1,
		uf_.fnp1+m-M6,	uf_.fn+m-M6,	uf_.fnm1+m-M6,
		uf_.fn  +m-M1,	uf_.fn+m-M2,	uf_.fn  +m-M4,	uf_.fn+m-M8,
		uf_.fn  +m+M1,	uf_.fn+m-M3,	uf_.fn  +m-M5,	uf_.fn+m-M9);
	  }

	/* Loop over the edge points in the mesh on the z = (zmin,zmax) and y = (ymin,ymax) boundary and
	 * update the fields using the first order absorbing boundary condition. To understand the
	 * following lines of codes, it is better to list the values of Li at the side.			*/
	uf_.af.ufB_ = &uf_.fE[0];
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  {
	    if ( rank_ == 0 )
	      {
		m = N1_ * i;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+3,	uf_.an+l+3,	uf_.anm1+l+3,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l+L4,	uf_.an+l+L4,	uf_.anm1+l+L4,
		    uf_.an  +l-L0,	uf_.an+l+L7,	uf_.an  +l-L3,	uf_.an+l+L10,
		    uf_.an  +l+L0,	uf_.an+l+L6,	uf_.an  +l+L2,	uf_.an+l+L8);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m+1,	uf_.fn+m+1,	uf_.fnm1+m+1,
		    uf_.fnp1+m+M1,	uf_.fn+m+M1,	uf_.fnm1+m+M1,
		    uf_.fnp1+m+M4,	uf_.fn+m+M4,	uf_.fnm1+m+M4,
		    uf_.fn  +m-M0,	uf_.fn+m+M7,	uf_.fn  +m-M3,	uf_.fn+m+M10,
		    uf_.fn  +m+M0,	uf_.fn+m+M6,	uf_.fn  +m+M2,	uf_.fn+m+M8);

		m = N1_ * i + uf_.N1m1;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-3,	uf_.an+l-3,	uf_.anm1+l-3,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l-L5,	uf_.an+l-L5,	uf_.anm1+l-L5,
		    uf_.an  +l-L0,	uf_.an+l-L6,	uf_.an  +l-L3,	uf_.an+l-L9,
		    uf_.an  +l+L0,	uf_.an+l-L7,	uf_.an  +l+L2,	uf_.an+l-L11);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m-1,	uf_.fn+m-1,	uf_.fnm1+m-1,
		    uf_.fnp1+m+M1,	uf_.fn+m+M1,	uf_.fnm1+m+M1,
		    uf_.fnp1+m-M5,	uf_.fn+m-M5,	uf_.fnm1+m-M5,
		    uf_.fn  +m-M0,	uf_.fn+m-M6,	uf_.fn  +m-M3,	uf_.fn+m-M9,
		    uf_.fn  +m+M0,	uf_.fn+m-M7,	uf_.fn  +m+M2,	uf_.fn+m-M11);
	      }

	    if ( rank_ == size_ - 1 )
	      {
		m = N1_ * i + N1N0_ * uf_.npm1;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+3,	uf_.an+l+3,	uf_.anm1+l+3,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l+L5,	uf_.an+l+L5,	uf_.anm1+l+L5,
		    uf_.an  +l-L0,	uf_.an+l+L7,	uf_.an  +l-L2,	uf_.an+l+L11,
		    uf_.an  +l+L0,	uf_.an+l+L6,	uf_.an  +l+L3,	uf_.an+l+L9);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m+1,	uf_.fn+m+1,	uf_.fnm1+m+1,
		    uf_.fnp1+m-M1,	uf_.fn+m-M1,	uf_.fnm1+m-M1,
		    uf_.fnp1+m+M5,	uf_.fn+m+M5,	uf_.fnm1+m+M5,
		    uf_.fn  +m-M0,	uf_.fn+m+M7,	uf_.fn  +m-M2,	uf_.fn+m+M11,
		    uf_.fn  +m+M0,	uf_.fn+m+M6,	uf_.fn  +m+M3,	uf_.fn+m+M9);

		m = N1_ * i + N1N0_ * uf_.npm1 + uf_.N1m1;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-3,	uf_.an+l-3,	uf_.anm1+l-3,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l-L4,	uf_.an+l-L4,	uf_.anm1+l-L4,
		    uf_.an  +l-L0,	uf_.an+l-L6,	uf_.an  +l-L2,	uf_.an+l-L8,
		    uf_.an  +l+L0,	uf_.an+l-L7,	uf_.an  +l+L3,	uf_.an+l-L10);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m-1,	uf_.fn+m-1,	uf_.fnm1+m-1,
		    uf_.fnp1+m-M1,	uf_.fn+m-M1,	uf_.fnm1+m-M1,
		    uf_.fnp1+m-M4,	uf_.fn+m-M4,	uf_.fnm1+m-M4,
		    uf_.fn  +m-M0,	uf_.fn+m-M6,	uf_.fn  +m-M2,	uf_.fn+m-M8,
		    uf_.fn  +m+M0,	uf_.fn+m-M7,	uf_.fn  +m+M3,	uf_.fn+m-M10);
	      }
	  }

	/* Loop over the edge points in the mesh on the z = (zmin,zmax) and x = (xmin,xmax) boundary and
	 * update the fields using the first order absorbing boundary condition. To understand the
	 * following lines of codes, it is better to list the values of Li at the side.			*/
	uf_.af.ufB_ = &uf_.gE[0];
	for (unsigned j = 1; j < uf_.N1m1; j++)
	  {
	    if ( rank_ == 0 )
	      {
		m = j;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l+L0,	uf_.an+l+L0,	uf_.anm1+l+L0,
		    uf_.anp1+l+L2,	uf_.an+l+L2,	uf_.anm1+l+L2,
		    uf_.an  +l-3,	uf_.an+l-L5,	uf_.an  +l-L7,	uf_.an+l-L11,
		    uf_.an  +l+3,	uf_.an+l+L4,	uf_.an  +l+L6,	uf_.an+l+L8);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m+M1,	uf_.fn+m+M1,	uf_.fnm1+m+M1,
		    uf_.fnp1+m+M0,	uf_.fn+m+M0,	uf_.fnm1+m+M0,
		    uf_.fnp1+m+M2,	uf_.fn+m+M2,	uf_.fnm1+m+M2,
		    uf_.fn  +m-1,	uf_.fn+m-M5,	uf_.fn  +m-M7,	uf_.fn+m-M11,
		    uf_.fn  +m+1,	uf_.fn+m+M4,	uf_.fn  +m+M6,	uf_.fn+m+M8);

		m = N1N0_ - N1_ + j;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l-L0,	uf_.an+l-L0,	uf_.anm1+l-L0,
		    uf_.anp1+l-L3,	uf_.an+l-L3,	uf_.anm1+l-L3,
		    uf_.an  +l-3,	uf_.an+l-L5,	uf_.an  +l-L6,	uf_.an+l-L9,
		    uf_.an  +l+3,	uf_.an+l+L4,	uf_.an  +l+L7,	uf_.an+l+L10);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m+M1,	uf_.fn+m+M1,	uf_.fnm1+m+M1,
		    uf_.fnp1+m-M0,	uf_.fn+m-M0,	uf_.fnm1+m-M0,
		    uf_.fnp1+m-M3,	uf_.fn+m-M3,	uf_.fnm1+m-M3,
		    uf_.fn  +m-1,	uf_.fn+m-M5,	uf_.fn  +m-M6,	uf_.fn+m-M9,
		    uf_.fn  +m+1,	uf_.fn+m+M4,	uf_.fn  +m+M7,	uf_.fn+m+M10);
	      }

	    if ( rank_ == size_ - 1 )
	      {
		m = N1N0_ * uf_.npm1 + j;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l+L0,	uf_.an+l+L0,	uf_.anm1+l+L0,
		    uf_.anp1+l+L3,	uf_.an+l+L3,	uf_.anm1+l+L3,
		    uf_.an  +l-3,	uf_.an+l-L4,	uf_.an  +l-L7,	uf_.an+l-L10,
		    uf_.an  +l+3,	uf_.an+l+L5,	uf_.an  +l+L6,	uf_.an+l+L9);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m-M1,	uf_.fn+m-M1,	uf_.fnm1+m-M1,
		    uf_.fnp1+m+M0,	uf_.fn+m+M0,	uf_.fnm1+m+M0,
		    uf_.fnp1+m+M3,	uf_.fn+m+M3,	uf_.fnm1+m+M3,
		    uf_.fn  +m-1,	uf_.fn+m-M4,	uf_.fn  +m-M7,	uf_.fn+m-M10,
		    uf_.fn  +m+1,	uf_.fn+m+M5,	uf_.fn  +m+M6,	uf_.fn+m+M9);

		m = N1N0_ * uf_.npm1 + N1_ * uf_.N0m1 + j;
		l = 3 * m;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,	uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l-L0,	uf_.an+l-L0,	uf_.anm1+l-L0,
		    uf_.anp1+l-L2,	uf_.an+l-L2,	uf_.anm1+l-L2,
		    uf_.an  +l-3,	uf_.an+l-L4,	uf_.an  +l-L6,	uf_.an+l-L8,
		    uf_.an  +l+3,	uf_.an+l+L5,	uf_.an  +l+L7,	uf_.an+l+L11);
		uf_.af.advanceEdgeS(
		    uf_.fnp1+m,	uf_.fn+m,	uf_.fnm1+m,
		    uf_.fnp1+m-M1,	uf_.fn+m-M1,	uf_.fnm1+m-M1,
		    uf_.fnp1+m-M0,	uf_.fn+m-M0,	uf_.fnm1+m-M0,
		    uf_.fnp1+m-M2,	uf_.fn+m-M2,	uf_.fnm1+m-M2,
		    uf_.fn  +m-1,	uf_.fn+m-M4,	uf_.fn  +m-M6,	uf_.fn+m-M8,
		    uf_.fn  +m+1,	uf_.fn+m+M5,	uf_.fn  +m+M7,	uf_.fn+m+M11);
	      }
	  }

	/* Now update the fields of the eight corners in the computational domain.			*/
	uf_.af.ufB_ = &uf_.hC[0];
	if ( rank_ == 0 )
	  {
	    m = 0;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m+N1_+1),		uf_.an+3*(m+N1_+1),		uf_.anm1+3*(m+N1_+1),
		uf_.anp1+3*(m+N1N0_+N1_),	uf_.an+3*(m+N1N0_+N1_),		uf_.anm1+3*(m+N1N0_+N1_),
		uf_.anp1+3*(m+N1N0_+1),	uf_.an+3*(m+N1N0_+1),		uf_.anm1+3*(m+N1N0_+1),
		uf_.anp1+3*(m+N1N0_+N1_+1),	uf_.an+3*(m+N1N0_+N1_+1),	uf_.anm1+3*(m+N1N0_+N1_+1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m+N1_),		uf_.fn+(m+N1_),			uf_.fnm1+(m+N1_),
		uf_.fnp1+(m+1),		uf_.fn+(m+1),			uf_.fnm1+(m+1),
		uf_.fnp1+(m+N1N0_),		uf_.fn+(m+N1N0_),		uf_.fnm1+(m+N1N0_),
		uf_.fnp1+(m+N1_+1),		uf_.fn+(m+N1_+1),		uf_.fnm1+(m+N1_+1),
		uf_.fnp1+(m+N1N0_+N1_),	uf_.fn+(m+N1N0_+N1_),		uf_.fnm1+(m+N1N0_+N1_),
		uf_.fnp1+(m+N1N0_+1),	uf_.fn+(m+N1N0_+1),		uf_.fnm1+(m+N1N0_+1),
		uf_.fnp1+(m+N1N0_+N1_+1),	uf_.fn+(m+N1N0_+N1_+1),		uf_.fnm1+(m+N1N0_+N1_+1));

	    m = N1N0_ - N1_;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m-N1_+1),		uf_.an+3*(m-N1_+1),		uf_.anm1+3*(m-N1_+1),
		uf_.anp1+3*(m+N1N0_-N1_),	uf_.an+3*(m+N1N0_-N1_),		uf_.anm1+3*(m+N1N0_-N1_),
		uf_.anp1+3*(m+N1N0_+1),	uf_.an+3*(m+N1N0_+1),		uf_.anm1+3*(m+N1N0_+1),
		uf_.anp1+3*(m+N1N0_-N1_+1),	uf_.an+3*(m+N1N0_-N1_+1),	uf_.anm1+3*(m+N1N0_-N1_+1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m-N1_),		uf_.fn+(m-N1_),			uf_.fnm1+(m-N1_),
		uf_.fnp1+(m+1),		uf_.fn+(m+1),			uf_.fnm1+(m+1),
		uf_.fnp1+(m+N1N0_),		uf_.fn+(m+N1N0_),		uf_.fnm1+(m+N1N0_),
		uf_.fnp1+(m-N1_+1),		uf_.fn+(m-N1_+1),		uf_.fnm1+(m-N1_+1),
		uf_.fnp1+(m+N1N0_-N1_),	uf_.fn+(m+N1N0_-N1_),		uf_.fnm1+(m+N1N0_-N1_),
		uf_.fnp1+(m+N1N0_+1),	uf_.fn+(m+N1N0_+1),		uf_.fnm1+(m+N1N0_+1),
		uf_.fnp1+(m+N1N0_-N1_+1),	uf_.fn+(m+N1N0_-N1_+1),		uf_.fnm1+(m+N1N0_-N1_+1));

	    m = uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m+N1_-1),		uf_.an+3*(m+N1_-1),		uf_.anm1+3*(m+N1_-1),
		uf_.anp1+3*(m+N1N0_+N1_),	uf_.an+3*(m+N1N0_+N1_),		uf_.anm1+3*(m+N1N0_+N1_),
		uf_.anp1+3*(m+N1N0_-1),	uf_.an+3*(m+N1N0_-1),		uf_.anm1+3*(m+N1N0_-1),
		uf_.anp1+3*(m+N1N0_+N1_-1),	uf_.an+3*(m+N1N0_+N1_-1),	uf_.anm1+3*(m+N1N0_+N1_-1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m+N1_),		uf_.fn+(m+N1_),			uf_.fnm1+(m+N1_),
		uf_.fnp1+(m-1),		uf_.fn+(m-1),			uf_.fnm1+(m-1),
		uf_.fnp1+(m+N1N0_),		uf_.fn+(m+N1N0_),		uf_.fnm1+(m+N1N0_),
		uf_.fnp1+(m+N1_-1),		uf_.fn+(m+N1_-1),		uf_.fnm1+(m+N1_-1),
		uf_.fnp1+(m+N1N0_+N1_),	uf_.fn+(m+N1N0_+N1_),		uf_.fnm1+(m+N1N0_+N1_),
		uf_.fnp1+(m+N1N0_-1),	uf_.fn+(m+N1N0_-1),		uf_.fnm1+(m+N1N0_-1),
		uf_.fnp1+(m+N1N0_+N1_-1),	uf_.fn+(m+N1N0_+N1_-1),		uf_.fnm1+(m+N1N0_+N1_-1));

	    m = N1N0_ - N1_ + uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m-N1_-1),		uf_.an+3*(m-N1_-1),		uf_.anm1+3*(m-N1_-1),
		uf_.anp1+3*(m+N1N0_-N1_),	uf_.an+3*(m+N1N0_-N1_),		uf_.anm1+3*(m+N1N0_-N1_),
		uf_.anp1+3*(m+N1N0_-1),	uf_.an+3*(m+N1N0_-1),		uf_.anm1+3*(m+N1N0_-1),
		uf_.anp1+3*(m+N1N0_-N1_-1),	uf_.an+3*(m+N1N0_-N1_-1),	uf_.anm1+3*(m+N1N0_-N1_-1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m-N1_),		uf_.fn+(m-N1_),			uf_.fnm1+(m-N1_),
		uf_.fnp1+(m-1),		uf_.fn+(m-1),			uf_.fnm1+(m-1),
		uf_.fnp1+(m+N1N0_),		uf_.fn+(m+N1N0_),		uf_.fnm1+(m+N1N0_),
		uf_.fnp1+(m-N1_-1),		uf_.fn+(m-N1_-1),		uf_.fnm1+(m-N1_-1),
		uf_.fnp1+(m+N1N0_-N1_),	uf_.fn+(m+N1N0_-N1_),		uf_.fnm1+(m+N1N0_-N1_),
		uf_.fnp1+(m+N1N0_-1),	uf_.fn+(m+N1N0_-1),		uf_.fnm1+(m+N1N0_-1),
		uf_.fnp1+(m+N1N0_-N1_-1),	uf_.fn+(m+N1N0_-N1_-1),		uf_.fnm1+(m+N1N0_-N1_-1));
	  }

	if ( rank_ == size_ - 1 )
	  {
	    m = N1N0_ * uf_.npm1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m+N1_+1),		uf_.an+3*(m+N1_+1),		uf_.anm1+3*(m+N1_+1),
		uf_.anp1+3*(m-N1N0_+N1_),	uf_.an+3*(m-N1N0_+N1_),		uf_.anm1+3*(m-N1N0_+N1_),
		uf_.anp1+3*(m-N1N0_+1),	uf_.an+3*(m-N1N0_+1),		uf_.anm1+3*(m-N1N0_+1),
		uf_.anp1+3*(m-N1N0_+N1_+1),	uf_.an+3*(m-N1N0_+N1_+1),	uf_.anm1+3*(m-N1N0_+N1_+1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m+N1_),		uf_.fn+(m+N1_),			uf_.fnm1+(m+N1_),
		uf_.fnp1+(m+1),		uf_.fn+(m+1),			uf_.fnm1+(m+1),
		uf_.fnp1+(m-N1N0_),		uf_.fn+(m-N1N0_),		uf_.fnm1+(m-N1N0_),
		uf_.fnp1+(m+N1_+1),		uf_.fn+(m+N1_+1),		uf_.fnm1+(m+N1_+1),
		uf_.fnp1+(m-N1N0_+N1_),	uf_.fn+(m-N1N0_+N1_),		uf_.fnm1+(m-N1N0_+N1_),
		uf_.fnp1+(m-N1N0_+1),	uf_.fn+(m-N1N0_+1),		uf_.fnm1+(m-N1N0_+1),
		uf_.fnp1+(m-N1N0_+N1_+1),	uf_.fn+(m-N1N0_+N1_+1),		uf_.fnm1+(m-N1N0_+N1_+1));

	    m = N1N0_ * uf_.npm1 + N1N0_ - N1_;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m-N1_+1),		uf_.an+3*(m-N1_+1),		uf_.anm1+3*(m-N1_+1),
		uf_.anp1+3*(m-N1N0_-N1_),	uf_.an+3*(m-N1N0_-N1_),		uf_.anm1+3*(m-N1N0_-N1_),
		uf_.anp1+3*(m-N1N0_+1),	uf_.an+3*(m-N1N0_+1),		uf_.anm1+3*(m-N1N0_+1),
		uf_.anp1+3*(m-N1N0_-N1_+1),	uf_.an+3*(m-N1N0_-N1_+1),	uf_.anm1+3*(m-N1N0_-N1_+1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m-N1_),		uf_.fn+(m-N1_),			uf_.fnm1+(m-N1_),
		uf_.fnp1+(m+1),		uf_.fn+(m+1),			uf_.fnm1+(m+1),
		uf_.fnp1+(m-N1N0_),		uf_.fn+(m-N1N0_),		uf_.fnm1+(m-N1N0_),
		uf_.fnp1+(m-N1_+1),		uf_.fn+(m-N1_+1),		uf_.fnm1+(m-N1_+1),
		uf_.fnp1+(m-N1N0_-N1_),	uf_.fn+(m-N1N0_-N1_),		uf_.fnm1+(m-N1N0_-N1_),
		uf_.fnp1+(m-N1N0_+1),	uf_.fn+(m-N1N0_+1),		uf_.fnm1+(m-N1N0_+1),
		uf_.fnp1+(m-N1N0_-N1_+1),	uf_.fn+(m-N1N0_-N1_+1),		uf_.fnm1+(m-N1N0_-N1_+1));

	    m = N1N0_ * uf_.npm1 + uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m+N1_-1),		uf_.an+3*(m+N1_-1),		uf_.anm1+3*(m+N1_-1),
		uf_.anp1+3*(m-N1N0_+N1_),	uf_.an+3*(m-N1N0_+N1_),		uf_.anm1+3*(m-N1N0_+N1_),
		uf_.anp1+3*(m-N1N0_-1),	uf_.an+3*(m-N1N0_-1),		uf_.anm1+3*(m-N1N0_-1),
		uf_.anp1+3*(m-N1N0_+N1_-1),	uf_.an+3*(m-N1N0_+N1_-1),	uf_.anm1+3*(m-N1N0_+N1_-1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m+N1_),		uf_.fn+(m+N1_),			uf_.fnm1+(m+N1_),
		uf_.fnp1+(m-1),		uf_.fn+(m-1),			uf_.fnm1+(m-1),
		uf_.fnp1+(m-N1N0_),		uf_.fn+(m-N1N0_),		uf_.fnm1+(m-N1N0_),
		uf_.fnp1+(m+N1_-1),		uf_.fn+(m+N1_-1),		uf_.fnm1+(m+N1_-1),
		uf_.fnp1+(m-N1N0_+N1_),	uf_.fn+(m-N1N0_+N1_),		uf_.fnm1+(m-N1N0_+N1_),
		uf_.fnp1+(m-N1N0_-1),	uf_.fn+(m-N1N0_-1),		uf_.fnm1+(m-N1N0_-1),
		uf_.fnp1+(m-N1N0_+N1_-1),	uf_.fn+(m-N1N0_+N1_-1),		uf_.fnm1+(m-N1N0_+N1_-1));

	    m = N1N0_ * uf_.npm1 + N1N0_ - N1_ + uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m-N1_-1),		uf_.an+3*(m-N1_-1),		uf_.anm1+3*(m-N1_-1),
		uf_.anp1+3*(m-N1N0_-N1_),	uf_.an+3*(m-N1N0_-N1_),		uf_.anm1+3*(m-N1N0_-N1_),
		uf_.anp1+3*(m-N1N0_-1),	uf_.an+3*(m-N1N0_-1),		uf_.anm1+3*(m-N1N0_-1),
		uf_.anp1+3*(m-N1N0_-N1_-1),	uf_.an+3*(m-N1N0_-N1_-1),	uf_.anm1+3*(m-N1N0_-N1_-1));
	    uf_.af.advanceCornerS(
		uf_.fnp1+m,			uf_.fn+m,			uf_.fnm1+m,
		uf_.fnp1+(m-N1_),		uf_.fn+(m-N1_),			uf_.fnm1+(m-N1_),
		uf_.fnp1+(m-1),		uf_.fn+(m-1),			uf_.fnm1+(m-1),
		uf_.fnp1+(m-N1N0_),		uf_.fn+(m-N1N0_),		uf_.fnm1+(m-N1N0_),
		uf_.fnp1+(m-N1_-1),		uf_.fn+(m-N1_-1),		uf_.fnm1+(m-N1_-1),
		uf_.fnp1+(m-N1N0_-N1_),	uf_.fn+(m-N1N0_-N1_),		uf_.fnm1+(m-N1N0_-N1_),
		uf_.fnp1+(m-N1N0_-1),	uf_.fn+(m-N1N0_-1),		uf_.fnm1+(m-N1N0_-1),
		uf_.fnp1+(m-N1N0_-N1_-1),	uf_.fn+(m-N1N0_-N1_-1),		uf_.fnm1+(m-N1N0_-N1_-1));
	  }
      }

    /* Communicate the calculated fields throughout the processors.					*/
    if (rank_ != size_ - 1)
      {
	MPI_Send(uf_.anp1+3*(np_-2)*N1N0_, 	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag1,MPI_COMM_WORLD);
	MPI_Send(uf_.fnp1+(np_-2)*N1N0_,	  N1N0_,MPI_DOUBLE,rank_+1,msgtag2,MPI_COMM_WORLD);
      }

    if (rank_ != 0)
      {
	MPI_Recv(uf_.anp1,			3*N1N0_,MPI_DOUBLE,rank_-1,msgtag1,MPI_COMM_WORLD,&status);
	MPI_Recv(uf_.fnp1,		  	  N1N0_,MPI_DOUBLE,rank_-1,msgtag2,MPI_COMM_WORLD,&status);
      }

    if (rank_ != 0)
      {
	MPI_Send(uf_.anp1+3*N1N0_,	 	3*N1N0_,MPI_DOUBLE,rank_-1,msgtag3,MPI_COMM_WORLD);
	MPI_Send(uf_.fnp1+N1N0_,		  N1N0_,MPI_DOUBLE,rank_-1,msgtag4,MPI_COMM_WORLD);
      }

    if (rank_ != size_ - 1)
      {
	MPI_Recv(uf_.anp1+3*(np_-1)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag3,MPI_COMM_WORLD,&status);
	MPI_Recv(uf_.fnp1+(np_-1)*N1N0_,	  N1N0_,MPI_DOUBLE,rank_+1,msgtag4,MPI_COMM_WORLD,&status);
      }

    /* Now that A and phi quantities are updated, calculate E and B at boundary grid points for later
     * acceleration and power measurement.								*/
    for (unsigned i = 1; i < uf_.N0m1; i++)
      for (unsigned j = 1; j < uf_.N1m1; j++)
	{
	  m = N1N0_ + N1_ * i + j;

	  /* Evaluate the field at this pixel.                                                        	*/
	  fieldEvaluate(m);

	  /* Set the boolean flag for this pixel and the one before it to true.                       	*/
	  pic_[m-N1N0_] = true;

	  /* For the left boundary (z=zmin) just set the fields equal to the next z-plane.		*/
	  if ( rank_ == 0 )
	    {
	      en_[m-N1N0_] = en_[m];
	      bn_[m-N1N0_] = bn_[m];
	    }

	  m = N1N0_ * ( np_ - 2 ) + N1_ * i + j;

	  /* Evaluate the field at this pixel.                                                        	*/
	  fieldEvaluate(m);

	  /* Set the boolean flag for this pixel and the one in front of it to true.                  	*/
	  pic_[m+N1N0_] = true;

	  /* For the right boundary (z=zmax) just set the fields equal to the previous z-plane.		*/
	  if ( rank_ == size_ - 1 )
	    {
	      en_[m+N1N0_] = en_[m];
	      bn_[m+N1N0_] = bn_[m];
	    }
	}

    /* Communicate the calculated fields throughout the processors.					*/
    if (rank_ != size_ - 1)
      {
	MPI_Send(uf_.en+3*(np_-2)*N1N0_, 	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag5,MPI_COMM_WORLD);
	MPI_Send(uf_.bn+3*(np_-2)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag6,MPI_COMM_WORLD);
      }

    if (rank_ != 0)
      {
	MPI_Recv(uf_.en,			3*N1N0_,MPI_DOUBLE,rank_-1,msgtag5,MPI_COMM_WORLD,&status);
	MPI_Recv(uf_.bn,		  	3*N1N0_,MPI_DOUBLE,rank_-1,msgtag6,MPI_COMM_WORLD,&status);
      }

    if (rank_ != 0)
      {
	MPI_Send(uf_.en+3*N1N0_,	 	3*N1N0_,MPI_DOUBLE,rank_-1,msgtag7,MPI_COMM_WORLD);
	MPI_Send(uf_.bn+3*N1N0_,		3*N1N0_,MPI_DOUBLE,rank_-1,msgtag8,MPI_COMM_WORLD);
      }

    if (rank_ != size_ - 1)
      {
	MPI_Recv(uf_.en+3*(np_-1)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag7,MPI_COMM_WORLD,&status);
	MPI_Recv(uf_.bn+3*(np_-1)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag8,MPI_COMM_WORLD,&status);
      }
  }

  /******************************************************************************************************
   * Evaluate the field of the m'th pixel from the potentials.
   ******************************************************************************************************/

  void FdTdSC::fieldEvaluate (long int m)
  {
    /* Calculate the electric field.                                                                  	*/
    en_[m].dv ( - uf_.dt, (*anp1_)[m] );
    en_[m].mdv( - uf_.dt, (*an_)  [m] );

    en_[m][0] -= ( *(uf_.fn+m+N1_  ) - *(uf_.fn+m-N1_  ) ) / uf_.dx2;
    en_[m][1] -= ( *(uf_.fn+m+1    ) - *(uf_.fn+m-1    ) ) / uf_.dy2;
    en_[m][2] -= ( *(uf_.fn+m+N1N0_) - *(uf_.fn+m-N1N0_) ) / uf_.dz2;

    /* Calculate the magnetic field.                                                                  	*/
    bn_[m][0] = 0.5 * (
	( *(uf_.an  +3*(m+1    )+2 ) - *(uf_.an  +3*(m-1    )+2  ) ) / uf_.dy2 -
	( *(uf_.an  +3*(m+N1N0_)+1 ) - *(uf_.an  +3*(m-N1N0_)+1  ) ) / uf_.dz2 +
	( *(uf_.anp1+3*(m+1    )+2 ) - *(uf_.anp1+3*(m-1    )+2  ) ) / uf_.dy2 -
	( *(uf_.anp1+3*(m+N1N0_)+1 ) - *(uf_.anp1+3*(m-N1N0_)+1  ) ) / uf_.dz2 );

    bn_[m][1] = 0.5 * (
	( *(uf_.an  +3*(m+N1N0_)   ) - *(uf_.an  +3*(m-N1N0_)    ) ) / uf_.dz2 -
	( *(uf_.an  +3*(m+N1_  )+2 ) - *(uf_.an  +3*(m-N1_  )+2  ) ) / uf_.dx2 +
	( *(uf_.anp1+3*(m+N1N0_)   ) - *(uf_.anp1+3*(m-N1N0_)    ) ) / uf_.dz2 -
	( *(uf_.anp1+3*(m+N1_  )+2 ) - *(uf_.anp1+3*(m-N1_  )+2  ) ) / uf_.dx2 );

    bn_[m][2] = 0.5 * (
	( *(uf_.an  +3*(m+N1_  )+1 ) - *(uf_.an  +3*(m-N1_  )+1  ) ) / uf_.dx2 -
	( *(uf_.an  +3*(m+1    )   ) - *(uf_.an  +3*(m-1    )    ) ) / uf_.dy2 +
	( *(uf_.anp1+3*(m+N1_  )+1 ) - *(uf_.anp1+3*(m-N1_  )+1  ) ) / uf_.dx2 -
	( *(uf_.anp1+3*(m+1    )   ) - *(uf_.anp1+3*(m-1    )    ) ) / uf_.dy2 );

    /* Set the boolean flag of this pixel to true.                                                    	*/
    pic_[m] = true;
  }

  /******************************************************************************************************
   * Sample the field and save it to the given file.
   ******************************************************************************************************/

  void FdTdSC::fieldSample ()
  {
    /* Do the field sampling if and only if sampling points are residing within this processor range.	*/
    if ( sf_.N > 0 )
      {
	( *(sf_.file) ).setf(std::ios::scientific);
	( *(sf_.file) ).precision(4);

	/* Write time into the first column.                                                       	*/
	*(sf_.file) << time_ * gamma_ << "\t";

	for (unsigned int n = 0; n < sf_.N; ++n)
	  {
	    /* Get the position of the sampling point.							*/
	    sf_.position 	= seed_.samplingPosition_[n];

	    /* Get the indices of the sampling point.							*/
	    sf_.dxr = modf( ( sf_.position[0] - xmin_ ) / mesh_.meshResolution_[0] , &sf_.c1);
	    sf_.i   = (int) sf_.c1;
	    sf_.dyr = modf( ( sf_.position[1] - ymin_ ) / mesh_.meshResolution_[1] , &sf_.c1);
	    sf_.j   = (int) sf_.c1;
	    sf_.dzr = modf( ( sf_.position[2] - zmin_ ) / mesh_.meshResolution_[2] , &sf_.c1);
	    sf_.k   = (int) sf_.c1;
	    sf_.m   = ( sf_.k - k0_ ) * N1N0_ + sf_.i * N1_ + sf_.j;

	    /* Calculate the fields to find the values at the sampling point.				*/
	    if (!pic_[sf_.m            ])     fieldEvaluate(sf_.m            );
	    if (!pic_[sf_.m+N1_        ])     fieldEvaluate(sf_.m+N1_        );
	    if (!pic_[sf_.m+1          ])     fieldEvaluate(sf_.m+1          );
	    if (!pic_[sf_.m+N1_+1      ])     fieldEvaluate(sf_.m+N1_+1      );
	    if (!pic_[sf_.m+N1N0_      ])     fieldEvaluate(sf_.m+N1N0_      );
	    if (!pic_[sf_.m+N1N0_+N1_  ])     fieldEvaluate(sf_.m+N1N0_+N1_  );
	    if (!pic_[sf_.m+N1N0_+1    ])     fieldEvaluate(sf_.m+N1N0_+1    );
	    if (!pic_[sf_.m+N1N0_+N1_+1])     fieldEvaluate(sf_.m+N1N0_+N1_+1);

	    /* Calculate and interpolate the electric field to find the value at the sampling point.	*/
	    sf_.et.mv ((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), en_[sf_.m]);
	    sf_.et.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), en_[sf_.m+N1_]);
	    sf_.et.pmv((1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr), en_[sf_.m+1]);
	    sf_.et.pmv(sf_.dxr         * sf_.dyr           * (1.0 - sf_.dzr), en_[sf_.m+N1_+1]);
	    sf_.et.pmv((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr,         en_[sf_.m+N1N0_]);
	    sf_.et.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * sf_.dzr,         en_[sf_.m+N1N0_+N1_]);
	    sf_.et.pmv((1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr,         en_[sf_.m+N1N0_+1]);
	    sf_.et.pmv(sf_.dxr         * sf_.dyr           * sf_.dzr,         en_[sf_.m+N1N0_+N1_+1]);

	    /* Calculate and interpolate the magnetic field to find its value at the sampling point.	*/
	    sf_.bt.mv ((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), bn_[sf_.m]);
	    sf_.bt.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), bn_[sf_.m+N1_]);
	    sf_.bt.pmv((1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr), bn_[sf_.m+1]);
	    sf_.bt.pmv(sf_.dxr         * sf_.dyr           * (1.0 - sf_.dzr), bn_[sf_.m+N1_+1]);
	    sf_.bt.pmv((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr,         bn_[sf_.m+N1N0_]);
	    sf_.bt.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * sf_.dzr,         bn_[sf_.m+N1N0_+N1_]);
	    sf_.bt.pmv((1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr,         bn_[sf_.m+N1N0_+1]);
	    sf_.bt.pmv(sf_.dxr         * sf_.dyr           * sf_.dzr,         bn_[sf_.m+N1N0_+N1_+1]);

	    sf_.at.mv ((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), (*an_)[sf_.m]);
	    sf_.at.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), (*an_)[sf_.m+N1_]);
	    sf_.at.pmv((1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr), (*an_)[sf_.m+1]);
	    sf_.at.pmv(sf_.dxr         * sf_.dyr           * (1.0 - sf_.dzr), (*an_)[sf_.m+N1_+1]);
	    sf_.at.pmv((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr,         (*an_)[sf_.m+N1N0_]);
	    sf_.at.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * sf_.dzr,         (*an_)[sf_.m+N1N0_+N1_]);
	    sf_.at.pmv((1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr,         (*an_)[sf_.m+N1N0_+1]);
	    sf_.at.pmv(sf_.dxr         * sf_.dyr           * sf_.dzr,         (*an_)[sf_.m+N1N0_+N1_+1]);

	    sf_.jt.mv ((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), jn_[sf_.m]);
	    sf_.jt.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), jn_[sf_.m+N1_]);
	    sf_.jt.pmv((1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr), jn_[sf_.m+1]);
	    sf_.jt.pmv(sf_.dxr         * sf_.dyr           * (1.0 - sf_.dzr), jn_[sf_.m+N1_+1]);
	    sf_.jt.pmv((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr,         jn_[sf_.m+N1N0_]);
	    sf_.jt.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * sf_.dzr,         jn_[sf_.m+N1N0_+N1_]);
	    sf_.jt.pmv((1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr,         jn_[sf_.m+N1N0_+1]);
	    sf_.jt.pmv(sf_.dxr         * sf_.dyr           * sf_.dzr,         jn_[sf_.m+N1N0_+N1_+1]);

	    sf_.f  =   (1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr) * (*fn_)[sf_.m];
	    sf_.f +=	  sf_.dxr        * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr) * (*fn_)[sf_.m+N1_];
	    sf_.f +=	 (1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr) * (*fn_)[sf_.m+1];
	    sf_.f +=	  sf_.dxr        * sf_.dyr           * (1.0 - sf_.dzr) * (*fn_)[sf_.m+N1_+1];
	    sf_.f +=   (1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr         * (*fn_)[sf_.m+N1N0_];
	    sf_.f +=	  sf_.dxr        * (1.0 - sf_.dyr)   * sf_.dzr         * (*fn_)[sf_.m+N1N0_+N1_];
	    sf_.f +=	 (1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr         * (*fn_)[sf_.m+N1N0_+1];
	    sf_.f +=	  sf_.dxr        * sf_.dyr           * sf_.dzr         * (*fn_)[sf_.m+N1N0_+N1_+1];

	    sf_.q  = 	 (1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr) * rn_[sf_.m];
	    sf_.q +=	  sf_.dxr        * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr) * rn_[sf_.m+N1_];
	    sf_.q +=   (1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr) * rn_[sf_.m+1];
	    sf_.q +=	  sf_.dxr        * sf_.dyr           * (1.0 - sf_.dzr) * rn_[sf_.m+N1_+1];
	    sf_.q +=   (1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr         * rn_[sf_.m+N1N0_];
	    sf_.q +=    sf_.dxr        * (1.0 - sf_.dyr)   * sf_.dzr         * rn_[sf_.m+N1N0_+N1_];
	    sf_.q +=   (1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr         * rn_[sf_.m+N1N0_+1];
	    sf_.q +=    sf_.dxr        * sf_.dyr           * sf_.dzr         * rn_[sf_.m+N1N0_+N1_+1];

	    /* Write the coordinates in the next column.						*/
	    *(sf_.file) << sf_.position[0] << "\t";
	    *(sf_.file) << sf_.position[1] << "\t";
	    *(sf_.file) << sf_.position[2] << "\t";

	    /* Write the fields in the next columns.							*/
	    for (unsigned int i = 0; i < seed_.samplingField_.size(); i++)
	      {
		if 		( seed_.samplingField_[i] == Ex )
		  *(sf_.file) << gamma_ * sf_.et[0] + c0_ * sqrt( pow(gamma_, 2) - 1 ) * sf_.bt[1] << "\t";
		else if 	( seed_.samplingField_[i] == Ey )
		  *(sf_.file) << gamma_ * sf_.et[1] - c0_ * sqrt( pow(gamma_, 2) - 1 ) * sf_.bt[0] << "\t";
		else if	( seed_.samplingField_[i] == Ez )
		  *(sf_.file) << sf_.et[2] << "\t";

		else if	( seed_.samplingField_[i] == Bx )
		  *(sf_.file) << gamma_ * sf_.bt[0] - sqrt( pow(gamma_, 2) - 1 ) / c0_ * sf_.et[1] << "\t";
		else if	( seed_.samplingField_[i] == By )
		  *(sf_.file) << gamma_ * sf_.bt[1] + sqrt( pow(gamma_, 2) - 1 ) / c0_ * sf_.et[0] << "\t";
		else if	( seed_.samplingField_[i] == Bz )
		  *(sf_.file) << sf_.bt[2] << "\t";

		else if	( seed_.samplingField_[i] == Ax )
		  *(sf_.file) << sf_.at[0] << "\t";
		else if	( seed_.samplingField_[i] == Ay )
		  *(sf_.file) << sf_.at[1] << "\t";
		else if	( seed_.samplingField_[i] == Az )
		  *(sf_.file) << sf_.at[2] << "\t";

		else if	( seed_.samplingField_[i] == Jx )
		  *(sf_.file) << sf_.jt[0] << "\t";
		else if	( seed_.samplingField_[i] == Jy )
		  *(sf_.file) << sf_.jt[1] << "\t";
		else if	( seed_.samplingField_[i] == Jz )
		  *(sf_.file) << sf_.jt[2] << "\t";

		else if	( seed_.samplingField_[i] == F  )
		  *(sf_.file) << sf_.f << "\t";
		else if	( seed_.samplingField_[i] == Q  )
		  *(sf_.file) << sf_.q << "\t";
	      }
	  }

	/* Take the file cursor to the next line.							*/
	*(sf_.file) << std::endl;
      }
  }

  /******************************************************************************************************
   * Visualize the field as vtk files on the whole domain and save them to the file with given name.
   ******************************************************************************************************/

  void FdTdSC::fieldVisualizeAllDomain (unsigned int ivtk)
  {
    long int			m;

    /* The old files if existing should be deleted.                                          		*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /* Calculate the field to be visualized in the vtk files.						*/
    for (int k = 0; k < np_; k++ )
      for (int j = 1; j < N1_-1; j++)
	for (int i = 1; i < N0_-1; i++)
	  {
	    m = k * N1_ * N0_ + i * N1_ + j;

	    if (!pic_[m]) fieldEvaluate(m);

	    for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	      {
		if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[m][l] = en_[m][0];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[m][l] = en_[m][1];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[m][l] = en_[m][2];
		else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[m][l] = bn_[m][0];
		else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[m][l] = bn_[m][1];
		else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[m][l] = bn_[m][2];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[m][l] = (*an_)[m][0];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[m][l] = (*an_)[m][1];
		else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[m][l] = (*an_)[m][2];
		else if 	( seed_.vtk_[ivtk].field_[l] == Jx )	vf_[ivtk].v[m][l] = jn_[m][0];
		else if 	( seed_.vtk_[ivtk].field_[l] == Jy )	vf_[ivtk].v[m][l] = jn_[m][1];
		else if 	( seed_.vtk_[ivtk].field_[l] == Jz )	vf_[ivtk].v[m][l] = jn_[m][2];
		else if 	( seed_.vtk_[ivtk].field_[l] == Q  ) 	vf_[ivtk].v[m][l] = rn_[m];
		else if 	( seed_.vtk_[ivtk].field_[l] == F  ) 	vf_[ivtk].v[m][l] = (*fn_)[m];
	      }
	  }

    /* Write the initial data for the vtk file.                                                       */
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;

    /* Insert the coordinates of the grid for the charge points.                                      */
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++)
      for (int j = 0; j < N1_; j++)
	for (int i = 0; i < N0_; i++)
	  {
	    m = k * N1_ * N0_ + i * N1_ + j;
	    *vf_[ivtk].file << r_[m][0] << " " << r_[m][1] << " " << r_[m][2] 			<< std::endl;
	  }
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /* Insert each cell data into the vtk file.                                            		*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          	<< std::endl;

    /* Insert the point data based on the computed electric field.					*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    		<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" format=\"ascii\">"
	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ) ; k++ )
      for (int j = 0; j < N1_; j++)
	for (int i = 0; i < N0_; i++)
	  {
	    m = k * N1_ * N0_ + i * N1_ + j;
	    *vf_[ivtk].file << vf_[ivtk].v[m][0];
	    for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[m][l];
	    *vf_[ivtk].file << std::endl;
	  }
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    		<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /* Close the file.                                                                      		*/
    (*vf_[ivtk].file).close();

    /* Write the file connecting the parallel files.							*/

    if ( rank_ == 0)
      {
	/* Add the vtk suffix and the number of the vtk file to the file name.                     	*/
	vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-" + stringify(nTime_) + PTS_FILE_SUFFIX;
	unsigned int k0, np;

	/* It is assumed that all the files in the directory are deleted before running the code.    	*/
	vf_[ivtk].file = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

	*vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
	*vf_[ivtk].file << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" >"			<< std::endl;
	*vf_[ivtk].file << "<PStructuredGrid WholeExtent=\"0 " << N0_-1 << " 0 " << N1_-1 << " 0 "  <<
	    N2_-1 << "\" GhostLevel = \"0\" >"                                    			<< std::endl;

	/* Insert the coordinates of the grid for the charge cloud.                            		*/
	*vf_[ivtk].file << "<PPoints>"                                                          	<< std::endl;
	*vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"	<< std::endl;
	*vf_[ivtk].file << "</PPoints>"                                                             	<< std::endl;

	*vf_[ivtk].file << "<PPointData>"                                                           	<< std::endl;
	*vf_[ivtk].file << "<DataArray type=\"Float64\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" Name=\"field\" format=\"ascii\" />"
	    << std::endl;
	*vf_[ivtk].file << "</PPointData>"                                                          	<< std::endl;

	for (int i = 0; i < size_; ++i)
	  {
	    /* Evaluate the number of nodes in each processor.						*/
	    if ( size_ > 1 )
	      {
		if ( i == 0 )
		  {
		    np = N2_ / size_ + 1;
		    k0 = 0;
		  }
		else if ( i == size_ - 1 )
		  {
		    np = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 3;
		    k0 = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
		  }
		else
		  {
		    np = N2_ / size_ + 2;
		    k0 = i * ( N2_ / size_ ) - 1;
		  }
	      }
	    else
	      {
		np = N2_;
		k0 = 0;
	      }

	    vf_[ivtk].fileName = vf_[ivtk].name + "-p" + stringify(i) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
	    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_-1 << " 0 " << N1_-1 << " " <<
		k0 << " " << k0 + np - 2 + ( (i == size_ - 1) ? 1 : 0 ) << "\""
		<< " Source=\"" << vf_[ivtk].fileName << "\" />"                       			<< std::endl;
	  }
	*vf_[ivtk].file << "</PStructuredGrid>"                                                     	<< std::endl;
	*vf_[ivtk].file << "</VTKFile>"                                                             	<< std::endl;

	/* Close the file.										*/
	(*vf_[ivtk].file).close();
      }
  }

  /******************************************************************************************************
   * Visualize the field as vtk files in plane and save them to the file with the given name.
   ******************************************************************************************************/

  void FdTdSC::fieldVisualizeInPlane (unsigned int ivtk)
  {
    if	( seed_.vtk_[ivtk].plane_ == XNORMAL )  fieldVisualizeInPlaneXNormal(ivtk);
    else if   ( seed_.vtk_[ivtk].plane_ == YNORMAL )  fieldVisualizeInPlaneYNormal(ivtk);
    else if   ( seed_.vtk_[ivtk].plane_ == ZNORMAL )
      {
	if ( seed_.vtk_[ivtk].position_[2] < zp_[1] && seed_.vtk_[ivtk].position_[2] >= zp_[0] )
	  fieldVisualizeInPlaneZNormal(ivtk);
      }
  }

  /******************************************************************************************************
   * Visualize the field as vtk files in a plane normal to x axis and save them to the file with the
   * given name.
   ******************************************************************************************************/

  void FdTdSC::fieldVisualizeInPlaneXNormal (unsigned int ivtk)
  {
    unsigned int		i, n;
    long int			m;
    Double			dxr, c;

    /* The old files if existing should be deleted.                                          		*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /* Calculate the index of the cell at which the plane resides.					*/
    dxr = modf( ( seed_.vtk_[ivtk].position_[0] - xmin_ ) / mesh_.meshResolution_[0] , &c);
    i   = (int) c;

    /* Calculate the field to be visualized in the vtk files.						*/
    for (int k = 0; k < np_; k++ )
      for (int j = 1; j < N1_-1; j++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  n = k * N1_ + j;

	  if (!pic_[m]) 		fieldEvaluate(m);
	  if (!pic_[m + N1_]) 	fieldEvaluate(m + N1_);

	  for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	    {
	      if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[n][l] = en_[m][0] * ( 1.0 - dxr ) + en_[m + N1_][0] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[n][l] = en_[m][1] * ( 1.0 - dxr ) + en_[m + N1_][1] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[n][l] = en_[m][2] * ( 1.0 - dxr ) + en_[m + N1_][2] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[n][l] = bn_[m][0] * ( 1.0 - dxr ) + bn_[m + N1_][0] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[n][l] = bn_[m][1] * ( 1.0 - dxr ) + bn_[m + N1_][1] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[n][l] = bn_[m][2] * ( 1.0 - dxr ) + bn_[m + N1_][2] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[n][l] = (*an_)[m][0] * ( 1.0 - dxr ) + (*an_)[m + N1_][0] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[n][l] = (*an_)[m][1] * ( 1.0 - dxr ) + (*an_)[m + N1_][1] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[n][l] = (*an_)[m][2] * ( 1.0 - dxr ) + (*an_)[m + N1_][2] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jx )	vf_[ivtk].v[n][l] = jn_[m][0] * ( 1.0 - dxr ) + jn_[m + N1_][0] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jy )	vf_[ivtk].v[n][l] = jn_[m][1] * ( 1.0 - dxr ) + jn_[m + N1_][1] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jz )	vf_[ivtk].v[n][l] = jn_[m][2] * ( 1.0 - dxr ) + jn_[m + N1_][2] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Q  )	vf_[ivtk].v[n][l] = rn_[m]    * ( 1.0 - dxr ) + rn_[m + N1_]    * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == F  )	vf_[ivtk].v[n][l] = (*fn_)[m] * ( 1.0 - dxr ) + (*fn_)[m + N1_] * dxr;
	    }
	}

    /* Write the initial data for the vtk file.                                                       	*/
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0  0  0 " << N1_ - 1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\" 0 0 0 " << N1_-1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;

    /* Insert the coordinates of the grid for the charge points.                                      	*/
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++)
      for (int j = 0; j < N1_; j++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  *vf_[ivtk].file << r_[m][0] * ( 1.0 - dxr ) + r_[m + N1_][0] * dxr << " "
	      << r_[m][1] << " " << r_[m][2] 								<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /* Insert each cell data into the vtk file.                                            		*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          		<< std::endl;

    /* Insert the point data based on the computed electric field.					*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    		<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" format=\"ascii\">"
	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ) ; k++ )
      for (int j = 0; j < N1_; j++)
	{
	  n = k * N1_ + j;
	  *vf_[ivtk].file << vf_[ivtk].v[n][0];
	  for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[n][l];
	  *vf_[ivtk].file << std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    		<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /* Close the file.                                                                      		*/
    (*vf_[ivtk].file).close();

    /* Write the file connecting the parallel files.							*/

    if ( rank_ == 0)
      {
	/* Add the vtk suffix and the number of the vtk file to the file name.                     	*/
	vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-" + stringify(nTime_) + PTS_FILE_SUFFIX;
	unsigned int k0, np;

	/* It is assumed that all the files in the directory are deleted before running the code.    	*/
	vf_[ivtk].file = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

	*vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
	*vf_[ivtk].file << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" >"			<< std::endl;
	*vf_[ivtk].file << "<PStructuredGrid WholeExtent=\" 0 0 0 " << N1_-1 << " 0 "  <<
	    N2_-1 << "\" GhostLevel = \"0\" >"                                    			<< std::endl;

	/* Insert the coordinates of the grid for the charge cloud.                            		*/
	*vf_[ivtk].file << "<PPoints>"                                                          	<< std::endl;
	*vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"	<< std::endl;
	*vf_[ivtk].file << "</PPoints>"                                                             	<< std::endl;

	*vf_[ivtk].file << "<PPointData>"                                                           	<< std::endl;
	*vf_[ivtk].file << "<DataArray type=\"Float64\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" Name=\"field\" format=\"ascii\" />"
	    << std::endl;
	*vf_[ivtk].file << "</PPointData>"                                                          	<< std::endl;

	for (int i = 0; i < size_; ++i)
	  {
	    /* Evaluate the number of nodes in each processor.						*/
	    if ( size_ > 1 )
	      {
		if ( i == 0 )
		  {
		    np = N2_ / size_ + 1;
		    k0 = 0;
		  }
		else if ( i == size_ - 1 )
		  {
		    np = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 3;
		    k0 = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
		  }
		else
		  {
		    np = N2_ / size_ + 2;
		    k0 = i * ( N2_ / size_ ) - 1;
		  }
	      }
	    else
	      {
		np = N2_;
		k0 = 0;
	      }

	    vf_[ivtk].fileName = vf_[ivtk].name + "-p" + stringify(i) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
	    *vf_[ivtk].file << "<Piece Extent=\"0 0 0 " << N1_-1 << " " <<
		k0 << " " << k0 + np - 2 + ( (i == size_ - 1) ? 1 : 0 ) << "\""
		<< " Source=\"" << vf_[ivtk].fileName << "\" />"                       			<< std::endl;
	  }
	*vf_[ivtk].file << "</PStructuredGrid>"                                                     	<< std::endl;
	*vf_[ivtk].file << "</VTKFile>"                                                             	<< std::endl;

	/* Close the file.										*/
	(*vf_[ivtk].file).close();
      }
  }

  /******************************************************************************************************
   * Visualize the field as vtk files in a plane normal to y axis and save them to the file with the
   * given name.
   ******************************************************************************************************/

  void FdTdSC::fieldVisualizeInPlaneYNormal (unsigned int ivtk)
  {
    unsigned int		j, n;
    long int			m;
    Double			dyr, c;

    /* The old files if existing should be deleted.                                          		*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /* Calculate the index of the cell at which the plane resides.					*/
    dyr = modf( ( seed_.vtk_[ivtk].position_[1] - ymin_ ) / mesh_.meshResolution_[1] , &c);
    j   = (int) c;

    /* Calculate the field to be visualized in the vtk files.						*/
    for (int k = 0; k < np_; k++ )
      for (int i = 1; i < N0_-1; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  n = k * N0_ + i;

	  if (!pic_[m]) 	fieldEvaluate(m);
	  if (!pic_[m + 1]) 	fieldEvaluate(m + 1);

	  for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	    {
	      if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[n][l] = en_[m][0] * ( 1.0 - dyr ) + en_[m + 1][0] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[n][l] = en_[m][1] * ( 1.0 - dyr ) + en_[m + 1][1] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[n][l] = en_[m][2] * ( 1.0 - dyr ) + en_[m + 1][2] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[n][l] = bn_[m][0] * ( 1.0 - dyr ) + bn_[m + 1][0] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[n][l] = bn_[m][1] * ( 1.0 - dyr ) + bn_[m + 1][1] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[n][l] = bn_[m][2] * ( 1.0 - dyr ) + bn_[m + 1][2] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[n][l] = (*an_)[m][0] * ( 1.0 - dyr ) + (*an_)[m + 1][0] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[n][l] = (*an_)[m][1] * ( 1.0 - dyr ) + (*an_)[m + 1][1] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[n][l] = (*an_)[m][2] * ( 1.0 - dyr ) + (*an_)[m + 1][2] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jx )	vf_[ivtk].v[n][l] = jn_[m][0] * ( 1.0 - dyr ) + jn_[m + 1][0] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jy )	vf_[ivtk].v[n][l] = jn_[m][1] * ( 1.0 - dyr ) + jn_[m + 1][1] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jz )	vf_[ivtk].v[n][l] = jn_[m][2] * ( 1.0 - dyr ) + jn_[m + 1][2] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Q  )	vf_[ivtk].v[n][l] = rn_[m]    * ( 1.0 - dyr ) + rn_[m + 1]    * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == F  )	vf_[ivtk].v[n][l] = (*fn_)[m] * ( 1.0 - dyr ) + (*fn_)[m + 1] * dyr;
	    }
	}

    /* Write the initial data for the vtk file.                                                       */
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 0 " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 0 " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;

    /* Insert the coordinates of the grid for the charge points.                                      */
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++)
      for (int i = 0; i < N0_; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  *vf_[ivtk].file << r_[m][0] * ( 1.0 - dyr ) + r_[m + 1][0] * dyr << " "
	      << r_[m][1] << " " << r_[m][2] 			<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /* Insert each cell data into the vtk file.                                            		*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          		<< std::endl;

    /* Insert the point data based on the computed electric field.					*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    	<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" format=\"ascii\">"
	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ) ; k++ )
      for (int i = 0; i < N0_; i++)
	{
	  n = k * N0_ + i;
	  *vf_[ivtk].file << vf_[ivtk].v[n][0];
	  for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[n][l];
	  *vf_[ivtk].file << std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    		<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /* Close the file.                                                                      		*/
    (*vf_[ivtk].file).close();

    /* Write the file connecting the parallel files.							*/

    if ( rank_ == 0)
      {
	/* Add the vtk suffix and the number of the vtk file to the file name.                     	*/
	vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-" + stringify(nTime_) + PTS_FILE_SUFFIX;
	unsigned int k0, np;

	/* It is assumed that all the files in the directory are deleted before running the code.    	*/
	vf_[ivtk].file = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

	*vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
	*vf_[ivtk].file << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" >"			<< std::endl;
	*vf_[ivtk].file << "<PStructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 0 0 "  <<
	    N2_-1 << "\" GhostLevel = \"0\" >"                                    			<< std::endl;

	/* Insert the coordinates of the grid for the charge cloud.                            		*/
	*vf_[ivtk].file << "<PPoints>"                                                          	<< std::endl;
	*vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"	<< std::endl;
	*vf_[ivtk].file << "</PPoints>"                                                             	<< std::endl;

	*vf_[ivtk].file << "<PPointData>"                                                           	<< std::endl;
	*vf_[ivtk].file << "<DataArray type=\"Float64\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" Name=\"field\" format=\"ascii\" />"
	    << std::endl;
	*vf_[ivtk].file << "</PPointData>"                                                          	<< std::endl;

	for (int i = 0; i < size_; ++i)
	  {
	    /* Evaluate the number of nodes in each processor.						*/
	    if ( size_ > 1 )
	      {
		if ( i == 0 )
		  {
		    np = N2_ / size_ + 1;
		    k0 = 0;
		  }
		else if ( i == size_ - 1 )
		  {
		    np = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 3;
		    k0 = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
		  }
		else
		  {
		    np = N2_ / size_ + 2;
		    k0 = i * ( N2_ / size_ ) - 1;
		  }
	      }
	    else
	      {
		np = N2_;
		k0 = 0;
	      }

	    vf_[ivtk].fileName = vf_[ivtk].name + "-p" + stringify(i) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
	    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 0 " <<
		k0 << " " << k0 + np - 2 + ( (i == size_ - 1) ? 1 : 0 ) << "\""
		<< " Source=\"" << vf_[ivtk].fileName << "\" />"                       			<< std::endl;
	  }
	*vf_[ivtk].file << "</PStructuredGrid>"                                                     	<< std::endl;
	*vf_[ivtk].file << "</VTKFile>"                                                             	<< std::endl;

	/* Close the file.										*/
	(*vf_[ivtk].file).close();
      }
  }

  /******************************************************************************************************
   * Visualize the field as vtk files in a plane normal to z axis and save them to the file with the
   * given name.
   ******************************************************************************************************/

  void FdTdSC::fieldVisualizeInPlaneZNormal (unsigned int ivtk)
  {
    unsigned int		k, n;
    long int			m;
    Double			dzr, c;

    /* The old files if existing should be deleted.                                          		*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /* Calculate the index of the cell at which the plane resides.					*/
    dzr = modf( ( seed_.vtk_[ivtk].position_[2] - zmin_ ) / mesh_.meshResolution_[2] , &c);
    k   = (int) c - k0_;

    /* Calculate the field to be visualized in the vtk files.						*/
    for (int j = 1; j < N1_-1; j++)
      for (int i = 1; i < N0_-1; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  n = i * N1_ + j;

	  if (!pic_[m]) 		fieldEvaluate(m);
	  if (!pic_[m + N1N0_]) 	fieldEvaluate(m + N1N0_);

	  for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	    {
	      if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[n][l] = en_[m][0] * ( 1.0 - dzr ) + en_[m + N1N0_][0] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[n][l] = en_[m][1] * ( 1.0 - dzr ) + en_[m + N1N0_][1] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[n][l] = en_[m][2] * ( 1.0 - dzr ) + en_[m + N1N0_][2] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[n][l] = bn_[m][0] * ( 1.0 - dzr ) + bn_[m + N1N0_][0] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[n][l] = bn_[m][1] * ( 1.0 - dzr ) + bn_[m + N1N0_][1] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[n][l] = bn_[m][2] * ( 1.0 - dzr ) + bn_[m + N1N0_][2] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[n][l] = (*an_)[m][0] * ( 1.0 - dzr ) + (*an_)[m + N1N0_][0] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[n][l] = (*an_)[m][1] * ( 1.0 - dzr ) + (*an_)[m + N1N0_][1] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[n][l] = (*an_)[m][2] * ( 1.0 - dzr ) + (*an_)[m + N1N0_][2] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jx )	vf_[ivtk].v[n][l] = jn_[m][0] * ( 1.0 - dzr ) + jn_[m + N1N0_][0] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jy )	vf_[ivtk].v[n][l] = jn_[m][1] * ( 1.0 - dzr ) + jn_[m + N1N0_][1] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Jz )	vf_[ivtk].v[n][l] = jn_[m][2] * ( 1.0 - dzr ) + jn_[m + N1N0_][2] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Q  )	vf_[ivtk].v[n][l] = rn_[m]    * ( 1.0 - dzr ) + rn_[m + N1N0_]    * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == F  )	vf_[ivtk].v[n][l] = (*fn_)[m] * ( 1.0 - dzr ) + (*fn_)[m + N1N0_] * dzr;
	    }
	}

    /* Write the initial data for the vtk file.                                                       	*/
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	0 << " " << 0 << "\">"										<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	0 << " " << 0 << "\">"										<< std::endl;

    /* Insert the coordinates of the grid for the charge points.                                      	*/
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    for (int j = 0; j < N1_; j++)
      for (int i = 0; i < N0_; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  *vf_[ivtk].file << r_[m][0] * ( 1.0 - dzr ) + r_[m + N1N0_][0] * dzr << " "
	      << r_[m][1] << " " << r_[m][2] 								<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /* Insert each cell data into the vtk file.                                            		*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          		<< std::endl;

    /* Insert the point data based on the computed electric field.					*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    		<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size()
				  << "\" format=\"ascii\">"						<< std::endl;
    for (int j = 0; j < N1_; j++)
      for (int i = 0; i < N0_; i++)
	{
	  n = i * N1_ + j;
	  *vf_[ivtk].file << vf_[ivtk].v[n][0];
	  for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[n][l];
	  *vf_[ivtk].file 										<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    		<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /* Close the file.                                                                      		*/
    (*vf_[ivtk].file).close();
  }

  /******************************************************************************************************
   * Write the total profile of the field into the given file name.
   ******************************************************************************************************/

  void FdTdSC::fieldProfile ()
  {
    /* Declare the iterators for the loop over the points.						*/
    pf_.fileName = seed_.profileBasename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + TXT_FILE_SUFFIX;
    pf_.file = new std::ofstream(pf_.fileName.c_str(),std::ios::trunc);

    pf_.file->setf(std::ios::scientific);
    pf_.file->precision(4);

    /* Perform a loop over the points of the mesh and save the field data into a text file.		*/
    for (int i = 0; i < N0_; i++ )
      for (int j = 0; j < N1_; j++ )
	for (int k = ( (rank_ == 0) ? 0 : 1 ); k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++ )
	  {
	    pf_.m = k * N1_ * N0_ + i * N1_ + j;
	    *pf_.file << r_[pf_.m][0] << "\t" ;
	    *pf_.file << r_[pf_.m][1] << "\t" ;
	    *pf_.file << r_[pf_.m][2] << "\t" ;

	    for (unsigned l = 0; l < seed_.profileField_.size(); l++)
	      {
		if 		( seed_.profileField_[l] == Ex )
		  *pf_.file << en_[pf_.m][0] << "\t";
		else if	( seed_.profileField_[l] == Ey )
		  *pf_.file << en_[pf_.m][1] << "\t";
		else if	( seed_.profileField_[l] == Ez )
		  *pf_.file << en_[pf_.m][2] << "\t";

		else if	( seed_.profileField_[l] == Bx )
		  *pf_.file << bn_[pf_.m][0] << "\t";
		else if	( seed_.profileField_[l] == By )
		  *pf_.file << bn_[pf_.m][1] << "\t";
		else if	( seed_.profileField_[l] == Bz )
		  *pf_.file << bn_[pf_.m][2] << "\t";

		else if 	( seed_.profileField_[l] == Ax )
		  *pf_.file << (*an_)[pf_.m][0] << "\t";
		else if 	( seed_.profileField_[l] == Ay )
		  *pf_.file << (*an_)[pf_.m][1] << "\t";
		else if 	( seed_.profileField_[l] == Az )
		  *pf_.file << (*an_)[pf_.m][2] << "\t";

		else if 	( seed_.profileField_[l] == Jx )
		  *pf_.file << jn_[pf_.m][0] << "\t";
		else if 	( seed_.profileField_[l] == Jy )
		  *pf_.file << jn_[pf_.m][1] << "\t";
		else if 	( seed_.profileField_[l] == Jz )
		  *pf_.file << jn_[pf_.m][2] << "\t";

		else if 	( seed_.profileField_[l] == F  )
		  *pf_.file << (*fn_)[pf_.m] << "\t";
		else if 	( seed_.profileField_[l] == Q  )
		  *pf_.file << rn_[pf_.m] << "\t";
	      }

	    *pf_.file << std::endl;
	  }

    /* Close the file.											*/
    (*pf_.file).close();
  }
}
