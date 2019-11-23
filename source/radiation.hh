/********************************************************************************************************
 *  power.hh : Implementation of the functions for calculation of radiated power
 ********************************************************************************************************/

#ifndef POWER_HH_
#define POWER_HH_

namespace Darius
{

  /******************************************************************************************************
   * Initialize the data required for sampling and saving the radiation power at the given position.
   ******************************************************************************************************/

  void Solver::initializeRadiationPower()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the FEL radiation power data.") );
    rp_.clear(); rp_.resize(FEL_.size());

    /* Loop over the different FEL output parameters and initialize the power calculation if it is
     * activated.											*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* Initialize if and only if the sampling of the power is enabled.				*/
	if (!FEL_[jf].radiationPower_.sampling_) continue;

	/* Perform the lorentz boost for the sampling data.						*/
	for (unsigned int i = 0; i < FEL_[jf].radiationPower_.z_.size(); i++)
	  FEL_[jf].radiationPower_.z_[i] 	*= gamma_;
	FEL_[jf].radiationPower_.lineBegin_ 	*= gamma_;
	FEL_[jf].radiationPower_.lineEnd_   	*= gamma_;

	Double dl = fabs(FEL_[jf].radiationPower_.lineEnd_ - FEL_[jf].radiationPower_.lineBegin_) / FEL_[jf].radiationPower_.res_;

	/* If sampling type is plot over line initialize the positions according to the line begin and
	 * line end.                                                                              	*/
	if ( FEL_[jf].radiationPower_.samplingType_ == OVERLINE )
	  {
	    Double l = 0.0;
	    FieldVector<Double> position;
	    while ( fabs(l) < fabs(FEL_[jf].radiationPower_.lineEnd_ - FEL_[jf].radiationPower_.lineBegin_) )
	      {
		FEL_[jf].radiationPower_.z_.push_back( FEL_[jf].radiationPower_.lineBegin_ + l);
		l += dl;
	      }
	  }

	/* Set the number of sampling points in each processor.                                       	*/
	rp_[jf].N  = FEL_[jf].radiationPower_.z_.size();
	rp_[jf].Nz = 0;
	for (unsigned int i = 0; i < rp_[jf].N; i++)
	  if ( FEL_[jf].radiationPower_.z_[i] < zp_[1] && FEL_[jf].radiationPower_.z_[i] >= zp_[0] )
	    ++rp_[jf].Nz;

	/* Add the normalized wavelength sweep to the vector of wavelengths.				*/
	dl = ( FEL_[jf].radiationPower_.lambdaMax_ - FEL_[jf].radiationPower_.lambdaMin_ ) / FEL_[jf].radiationPower_.lambdaRes_;
	for (Double rw = FEL_[jf].radiationPower_.lambdaMin_; rw < FEL_[jf].radiationPower_.lambdaMax_; rw += dl)
	  FEL_[jf].radiationPower_.lambda_.push_back(rw);
	rp_[jf].Nl = FEL_[jf].radiationPower_.lambda_.size();

	/* Based on the number of points to calculate and the size of the threads, allocate memory for
	 * saving the powers.										*/
	rp_[jf].pL.resize(rp_[jf].Nl * rp_[jf].N, 0.0);
	rp_[jf].pG.resize(rp_[jf].Nl * rp_[jf].N, 0.0);

	rp_[jf].Nf = 0;

	/* Initialize the file streams to save the data.						*/
	rp_[jf].file.resize(rp_[jf].Nl);
	rp_[jf].w.resize(rp_[jf].Nl);
	for (unsigned int i = 0; i < rp_[jf].Nl; i++)
	  {
	    std::string baseFilename = "";
	    if (!(isabsolute(FEL_[jf].radiationPower_.basename_))) baseFilename = FEL_[jf].radiationPower_.directory_;
	    baseFilename += FEL_[jf].radiationPower_.basename_ + "-" + stringify(i) + TXT_FILE_SUFFIX;

	    /* If the directory of the baseFilename does not exist create this directory.		*/
	    createDirectory(baseFilename, rank_);

	    rp_[jf].file[i] = new std::ofstream(baseFilename.c_str(),std::ios::trunc);
	    ( *(rp_[jf].file[i]) ).setf(std::ios::scientific);
	    ( *(rp_[jf].file[i]) ).precision(15);
	    ( *(rp_[jf].file[i]) ).width(40);

	    /* Determine the number of time points needed to calculate the amplitude of each radiation
	     * harmonic. Here, we use the power inside three radiation cycles to calculate the
	     * Instantaneous power at the selected harmonic.						*/
	    Double dt = 3.0 * undulator_[0].lu_ / FEL_[jf].radiationPower_.lambda_[i] / ( gamma_ * c0_ );
	    rp_[jf].Nf = ( int( dt / mesh_.timeStep_ ) > rp_[jf].Nf ) ? int( dt / mesh_.timeStep_ ) : rp_[jf].Nf;

	    /* Calculate the angular frequency for each wavelength.					*/
	    rp_[jf].w[i] = 2 * PI / dt;
	  }

	/* Based on the obtained Nf, resize the vectors for saving the time domain data.		*/
	rp_[jf].fdt.resize(rp_[jf].Nf, std::vector<std::vector<Double> > (rp_[jf].Nz * N1_ * N0_, std::vector<Double> (4,0.0) ) );

	rp_[jf].dt  = mesh_.timeStep_;
	rp_[jf].dx  = mesh_.meshResolution_[0];
	rp_[jf].dy  = mesh_.meshResolution_[1];
	rp_[jf].dz  = mesh_.meshResolution_[2];

	rp_[jf].pc  = 2.0 * rp_[jf].dx * rp_[jf].dy / ( m0_ * rp_[jf].Nf * rp_[jf].Nf ) * pow(mesh_.lengthScale_,2) / pow(mesh_.timeScale_,3);
      }

    /* Loop over the different FEL output parameters and initialize the power calculation if it is
     * activated.											*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* Initialize if and only if the sampling of the power is enabled.				*/
	if (!FEL_[jf].vtk_.sampling_) continue;

	/* Return an error if the power visualization rhythm is still zero.				*/
	if ( FEL_[jf].vtk_.rhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The power visualization rhythm of the field is zero although power visualization is activated !!!") );
	    exit(1);
	  }

	/* Lorentz boost the power-visualization sampling rhythm to the electron rest frame.		*/
	FEL_[jf].vtk_.rhythm_		/= gamma_;

	/* Perform the Lorentz boost for the sampling data.						*/
	FEL_[jf].vtk_.z_ 		*= gamma_;

	/* Create the filename for saving the visualization data.					*/
	if (!(isabsolute(FEL_[jf].vtk_.basename_)))
	  FEL_[jf].vtk_.basename_ = FEL_[jf].vtk_.directory_ + FEL_[jf].vtk_.basename_;

	/* If the directory of the baseFilename does not exist create this directory.			*/
	createDirectory(FEL_[jf].vtk_.basename_, rank_);

	/* Set the number of sampling points in each processor.                                       	*/
	rp_[jf].N  = 1;
	rp_[jf].Nz = ( FEL_[jf].vtk_.z_ < zp_[1] && FEL_[jf].vtk_.z_ >= zp_[0] ) ? 1 : 0;
	rp_[jf].Nl = 1;

	/* Do not continue the loop if Nz is not equal to one.						*/
	if ( rp_[jf].Nz == 0 ) continue;

	/* Based on the number of points to calculate and the size of the threads, allocate memory for
	 * saving the powers.										*/
	rp_[jf].pL.resize(N1_*N0_, 0.0);
	rp_[jf].pG.clear();

	rp_[jf].Nf = 0;

	/* Initialize the file streams to save the data.						*/
	rp_[jf].file.resize(rp_[jf].Nz);
	rp_[jf].w.resize(rp_[jf].Nz);

	/* Determine the number of time points needed to calculate the amplitude of each radiation
	 * harmonic.											*/
	Double dt = undulator_[0].lu_ / FEL_[jf].vtk_.lambda_ / ( gamma_ * c0_ );
	rp_[jf].Nf = int( dt / mesh_.timeStep_ );

	/* Calculate the angular frequency for each wavelength.						*/
	rp_[jf].w[0] = 2 * PI / dt;

	/* Based on the obtained Nf, resize the vectors for saving the time domain data.		*/
	rp_[jf].fdt.resize(rp_[jf].Nf, std::vector<std::vector<Double> > (N1_ * N0_, std::vector<Double> (4,0.0) ) );

	rp_[jf].dt  = mesh_.timeStep_;
	rp_[jf].dx  = mesh_.meshResolution_[0];
	rp_[jf].dy  = mesh_.meshResolution_[1];
	rp_[jf].dz  = mesh_.meshResolution_[2];

	rp_[jf].pc  = 2.0 * rp_[jf].dx * rp_[jf].dy / ( m0_ * rp_[jf].Nf * rp_[jf].Nf ) * pow(mesh_.lengthScale_,2) / pow(mesh_.timeScale_,3);
      }

    printmessage(std::string(__FILE__), __LINE__, std::string(" The FEL radiation power data is initialized. :::") );
  }

  /******************************************************************************************************
   * Sample the radiation power at the given position and save it to the file.
   ******************************************************************************************************/

  void Solver::radiationPowerSample()
  {
    /* Declare the temporary parameters needed for calculating the radiated power.                    	*/
    unsigned int              	k, l, m, i, j, kz;
    long int                  	mi, ni;
    FieldVector<Double>       	et, bt;
    Complex                   	ew1, bw1, ew2, bw2;
    Double			ec, es, ea;

    /* Loop over the different FEL output parameters and calculate the radiation energy if the energy
     * calculation is activated.                                                                      	*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* Initialize if and only if the sampling of the power is enabled.                            	*/
	if (!FEL_[jf].radiationPower_.sampling_) continue;

	/* First reset all the previously calculated powers.                                          	*/
	for (k = 0; k < rp_[jf].N; ++k)
	  for (l = 0; l < rp_[jf].Nl; ++l)
	    rp_[jf].pL[k * rp_[jf].Nl + l] = 0.0;

	/* Set the index of the sampling point to zero.                                               	*/
	kz = 0;

	/* Loop over the sampling positions, transverse dicretizations, and frequency to calculate the
	 * radiated power at the specific point and frequency.                                        	*/

	/* k index loops over the sampling positions.							*/
	for (k = 0; k < rp_[jf].N; ++k)
	  {
	    /* Do not continue if this index is not supported by the processor.                       	*/
	    if ( !( FEL_[jf].radiationPower_.z_[k] < zp_[1] && FEL_[jf].radiationPower_.z_[k] >= zp_[0] ) ) continue;

	    /* Obtain the z index of the cell containing the point.                                   	*/
	    rp_[jf].dzr = modf( ( FEL_[jf].radiationPower_.z_[k] - zmin_ ) / mesh_.meshResolution_[2] , &rp_[jf].c);
	    rp_[jf].k   = (int) rp_[jf].c;

	    /* Loop over the transverse indices.                                                      	*/
	    for (i = 2; i < N0_ - 2; i += 1)
	      for (j = 2; j < N1_ - 2; j += 1)
		{
		  /* Get the index in the computation grid as well as the field storage grid.         	*/
		  mi = ( rp_[jf].k - k0_) * N1_* N0_ + i * N1_ + j;
		  ni = kz * N1_* N0_ + i * N1_ + j;

		  /* Evaluate the fields of the corresponding pixels.                                 	*/
		  if (!pic_[mi      ]) fieldEvaluate(mi      );
		  if (!pic_[mi+N1N0_]) fieldEvaluate(mi+N1N0_);

		  /* Calculate and interpolate the electric field to find the value at the bunch point.	*/
		  et[0] = ( 1.0 - rp_[jf].dzr ) * en_[mi][0] + rp_[jf].dzr * en_[mi+N1N0_][0];
		  et[1] = ( 1.0 - rp_[jf].dzr ) * en_[mi][1] + rp_[jf].dzr * en_[mi+N1N0_][1];

		  /* Calculate and interpolate the magnetic field to find its value at the bunch point.	*/
		  bt[0] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][0] + rp_[jf].dzr * bn_[mi+N1N0_][0];
		  bt[1] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][1] + rp_[jf].dzr * bn_[mi+N1N0_][1];

		  /* Transform the fields to the lab frame.                                           	*/
		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][0] = gamma_ * ( et[0] + c0_ * beta_ * bt[1] );
		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][1] = gamma_ * ( et[1] - c0_ * beta_ * bt[0] );

		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][2] = gamma_ * ( bt[0] - beta_ / c0_ * et[1] );
		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][3] = gamma_ * ( bt[1] + beta_ / c0_ * et[0] );

		  /* Add the contribution of this field to the radiation power.                       	*/

		  /* l index loops over the given wavelength of the power sampling.			*/
		  for ( l = 0; l < rp_[jf].Nl; l++)
		    {
		      ew1 = Complex (0.0, 0.0);
		      bw1 = ew1;
		      ew2 = ew1;
		      bw2 = ew1;

		      for ( m = 0; m < rp_[jf].Nf; m++)
			{
			  ea	 = rp_[jf].w[l] * m * mesh_.timeStep_;
			  ec   = cos(ea);
			  es   = sin(ea);
			  ew1 += rp_[jf].fdt[m][ni][0] * ( ec + I * es );
			  bw1 += rp_[jf].fdt[m][ni][3] * ( ec - I * es );
			  ew2 += rp_[jf].fdt[m][ni][1] * ( ec + I * es );
			  bw2 += rp_[jf].fdt[m][ni][2] * ( ec - I * es );
			}

		      /* Add the contribution to the power series.					*/
		      rp_[jf].pL[k * rp_[jf].Nl + l] += rp_[jf].pc * ( std::real( ew1 * bw1 ) - std::real( ew2 * bw2 ) );
		    }
		}

	    /* Add the iterator for the sampling point by one.                                       	*/
	    kz += 1;
	  }

	/* Add the data from each processor together at the root processor.                           	*/
	MPI_Allreduce(&rp_[jf].pL[0],&rp_[jf].pG[0],rp_[jf].N*rp_[jf].Nl,MPI_TYPE,MPI_SUM,MPI_COMM_WORLD);

	/* If the rank of the processor is equal to zero, i.e. root processor save the fields into the
	 * given file.                                                                                	*/
	for ( l = 0; l < rp_[jf].Nl; l++)
	  {
	    if ( rank_ == ( l % size_ ) )
	      {
		for (k = 0; k < rp_[jf].N; ++k)
		  *(rp_[jf].file[l]) << gamma_ * ( FEL_[jf].radiationPower_.z_[k] + beta_ * c0_ * ( timeBunch_ + dt_ ) ) << "\t" << rp_[jf].pG[k * rp_[jf].Nl + l] << "\t";
		*(rp_[jf].file[l]) << std::endl;
	      }
	  }
      }
  }

  /******************************************************************************************************
   * Visualize the radiation power at the given position and save it to the file.
   ******************************************************************************************************/

  void Solver::radiationPowerVisualize()
  {

    /* Declare the temporary parameters needed for calculating the radiated power.                    	*/
    unsigned int              	k, l, m, i, j, kz;
    long int                  	mi, ni;
    FieldVector<Double>       	et, bt;
    Complex                   	ew1, bw1, ew2, bw2;
    Double			ec, es, ea;

    /* Loop over the different FEL output parameters and visualize the radiation power if visualization
     * is activated.                                                                      		*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {

	/* Initialize if and only if the sampling of the power is enabled.                            	*/
	if ( FEL_[jf].vtk_.sampling_ && rp_[jf].Nz == 1 )
	  {

	    /* Calculate the index of the cell at which the plane resides.				*/
	    rp_[jf].dzr 	= modf( ( FEL_[jf].vtk_.z_ - zmin_ ) / mesh_.meshResolution_[2] , &rp_[jf].c);
	    rp_[jf].k   	= (int) rp_[jf].c - k0_;

	    /* Loop over the transverse indices.                                                      	*/
	    for (i = 1; i < N0_ - 1; i++)
	      for (j = 1; j < N1_ - 1; j++)
		{
		  /* Get the index in the computation grid as well as the field storage grid.         	*/
		  mi = rp_[jf].k * N1_* N0_ + i * N1_ + j;
		  ni = i * N1_ + j;

		  /* Evaluate the fields of the corresponding pixels.                                 	*/
		  if (!pic_[mi      ]) fieldEvaluate(mi      );
		  if (!pic_[mi+N1N0_]) fieldEvaluate(mi+N1N0_);

		  /* Calculate and interpolate the electric field to find the value at the bunch point.	*/
		  et[0] = ( 1.0 - rp_[jf].dzr ) * en_[mi][0] + rp_[jf].dzr * en_[mi+N1N0_][0];
		  et[1] = ( 1.0 - rp_[jf].dzr ) * en_[mi][1] + rp_[jf].dzr * en_[mi+N1N0_][1];

		  /* Calculate and interpolate the magnetic field to find its value at the bunch point.	*/
		  bt[0] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][0] + rp_[jf].dzr * bn_[mi+N1N0_][0];
		  bt[1] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][1] + rp_[jf].dzr * bn_[mi+N1N0_][1];

		  /* Transform the fields to the lab frame.                                           	*/
		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][0] = gamma_ * ( et[0] + c0_ * beta_ * bt[1] );
		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][1] = gamma_ * ( et[1] - c0_ * beta_ * bt[0] );

		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][2] = gamma_ * ( bt[0] - beta_ / c0_ * et[1] );
		  rp_[jf].fdt[nTime_ % rp_[jf].Nf][ni][3] = gamma_ * ( bt[1] + beta_ / c0_ * et[0] );

		  /* Add the contribution of this field to the radiation power.                       	*/
		  ew1 = Complex (0.0, 0.0);
		  bw1 = ew1;
		  ew2 = ew1;
		  bw2 = ew1;

		  for ( m = 0; m < rp_[jf].Nf; m++)
		    {
		      ea   = rp_[jf].w[0] * m * mesh_.timeStep_;
		      ec   = cos(ea);
		      es   = sin(ea);
		      ew1 += rp_[jf].fdt[m][ni][0] * ( ec + I * es );
		      bw1 += rp_[jf].fdt[m][ni][3] * ( ec - I * es );
		      ew2 += rp_[jf].fdt[m][ni][1] * ( ec + I * es );
		      bw2 += rp_[jf].fdt[m][ni][2] * ( ec - I * es );
		    }

		  rp_[jf].pL[ni] = rp_[jf].pc * ( std::real( ew1 * bw1 ) - std::real( ew2 * bw2 ) );
		}

	    if ( fmod(time_, FEL_[jf].vtk_.rhythm_) < mesh_.timeStep_ )
	      {

		/* The old files if existing should be deleted.						*/
		std::string baseFilename = FEL_[jf].vtk_.basename_ + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
		rp_[jf].file[0] = new std::ofstream(baseFilename.c_str(),std::ios::trunc);

		rp_[jf].file[0]->setf(std::ios::scientific);
		rp_[jf].file[0]->precision(4);

		/* Write the initial data for the vtk file.						*/
		*rp_[jf].file[0] << "<?xml version=\"1.0\"?>"						<< std::endl;
		*rp_[jf].file[0] << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
		    "compressor=\"vtkZLibDataCompressor\">" 						<< std::endl;
		*rp_[jf].file[0] << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
		    0 << " " << 0 << "\">"								<< std::endl;
		*rp_[jf].file[0] << "<Piece Extent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
		    0 << " " << 0 << "\">"								<< std::endl;

		/* Insert the coordinates of the grid for the charge points.				*/
		*rp_[jf].file[0] << "<Points>"                                                        << std::endl;
		*rp_[jf].file[0] << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<< std::endl;
		for (j = 0; j < N1_; j++)
		  for (i = 0; i < N0_; i++)
		    {
		      m = rp_[jf].k * N1_ * N0_ + i * N1_ + j;
		      *rp_[jf].file[0] << r_[m][0] * ( 1.0 - rp_[jf].dzr ) + r_[m + N1N0_][0] * rp_[jf].dzr << " "
			  << r_[m][1] << " " << r_[m][2] 						<< std::endl;
		    }
		*rp_[jf].file[0] << "</DataArray>"                                                    << std::endl;
		*rp_[jf].file[0] << "</Points>"                                                       << std::endl;

		/* Insert each cell data into the vtk file.                                         	*/
		*rp_[jf].file[0] << "<CellData>"                                                       << std::endl;
		*rp_[jf].file[0] << "</CellData>"                                                      << std::endl;

		/* Insert the point data based on the computed electric field.				*/
		*rp_[jf].file[0] << "<PointData Vectors = \"power\">"                                 << std::endl;
		*rp_[jf].file[0] << "<DataArray type=\"Float64\" Name=\"power\" NumberOfComponents=\"" << 1 << "\" format=\"ascii\">"
		    << std::endl;
		for (j = 0; j < N1_; j++)
		  for (i = 0; i < N0_; i++)
		    *rp_[jf].file[0] << rp_[jf].pL[i * N1_ + j] 					<< std::endl;

		*rp_[jf].file[0] << "</DataArray>"                                                   	<< std::endl;
		*rp_[jf].file[0] << "</PointData>"                                                    << std::endl;
		*rp_[jf].file[0] << "</Piece>"                                                        << std::endl;
		*rp_[jf].file[0] << "</StructuredGrid>"                                               << std::endl;
		*rp_[jf].file[0] << "</VTKFile>"                                                      << std::endl;

		/* Close the file.                                                                    */
		(*rp_[jf].file[0]).close();
	      }
	  }
      }
  }

}       /* End of namespace Darius.                                                                    */
#endif
