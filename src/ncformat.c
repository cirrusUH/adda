/* File: ncformat.c
 * $Date::                            $
 * Descr:
 *
 *
 * Copyright (C) 2016 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef NETCDF4
  #error This should have been excluded from compiling, as NETCDF4 is not #defined.
#endif

#include "comm.h"
#include "io.h"
#include <stdlib.h>
#include <stdint.h>

#include "ncformat.h"

int g_ncid; // netcdf file handle
char g_ncmsgbuf[MAX_FNAME]; // buffer for filname (simplifies error msg handling)

static int g_geom_id; // shared in ncInitDipFile,ncReadDipFile

extern const int jagged;

// commandline parameters
extern bool netcdf4_addCFheader;

int ncFOpenErr(const char * restrict fname,int ncmode,ERR_LOC_DECL)
{
	int status;

#ifdef PARALLEL
	int mype,npe;

	MPI_Comm_rank (MPI_COMM_WORLD,&mype);
	MPI_Comm_size (MPI_COMM_WORLD,&npe);

	status=nc_open_par(fname,(ncmode|NC_NETCDF4|NC_MPIIO),MPI_COMM_WORLD,MPI_INFO_NULL,&g_ncid);
	// printf("ncFOpenErr: line=%d, status=%d, rank=%d, size=%d, err='%s'\n",__LINE__,status, mype, npe, nc_strerror(status) );

#else // SEQ

	status=nc_open(fname,ncmode|NC_NETCDF4,&g_ncid);
	// printf("ncFOpenErr: line=%d, status=%d,  err='%s'\n",__LINE__,status,  nc_strerror(status) );
#endif

	if (status!=NC_NOERR) {
		// printf("status=%d\n", status); // NC_ENOTNC == -51
		if (status==NC_ENOTNC) {
			LogWarning(EC_WARN,ERR_LOC_CALL,"'%s' is HDF5 but not a netCDF4 file.",fname);
			return status;
		 }

		 LogError(ERR_LOC_CALL,"Failed to open file '%s', ncerror='%s'",fname,nc_strerror(status));
	 }
	 LogWarning(EC_WARN,ONE_POS,"Opened NetCDF file '%s'",fname);
	 return g_ncid;
}

int ncFCreateErr(const char * restrict fname,int ncmode,ERR_LOC_DECL)
{
	int status;

#ifdef PARALLEL
	int mype,npe;

	MPI_Comm_rank (MPI_COMM_WORLD,&mype);
	MPI_Comm_size (MPI_COMM_WORLD,&npe);

	int mode=NC_MPIIO; // or NC_MPIPOSIX which is not preferred
	status=nc_create_par(fname,(ncmode|NC_NETCDF4|mode),MPI_COMM_WORLD,MPI_INFO_NULL,&g_ncid);

	// printf("ncFCreateErr: line=%d, status=%d, rank=%d, size=%d, err='%s'\n",__LINE__,status, mype, npe, nc_strerror(status) );
	if (status!=NC_NOERR) {
	  if (status==NC_ENOPAR) {
	     LogError(ERR_LOC_CALL,"Failed to create file '%s', ncerror='%s'. "\
		      "Check, using e.g. ldd ./adda_mpi, if the runtime linker finds the correct libnetcdf.so. "\
		      "Did you forget to set LD_LIBRARY_PATH?"\
		      ,fname,nc_strerror(status));
	  }
	  else {
		 LogError(ERR_LOC_CALL,"Failed to create file '%s', ncerror='%s'",fname,nc_strerror(status));
	  }
	}
#else // SEQ

	status=nc_create(fname,ncmode|NC_NETCDF4,&g_ncid);
	//printf("ncFCreateErr: line=%d, status=%d,  err='%s'\n",__LINE__,status,  nc_strerror(status) );
	if (status!=NC_NOERR) {
	     // printf("status=%d\n", status); // NC_ENOTNC == -51
		 LogError(ERR_LOC_CALL,"Failed to create file '%s', ncerror='%s'",fname,nc_strerror(status));
	 }

#endif

	LogWarning(EC_WARN,ONE_POS,"Created NetCDF file '%s'",ncGetFilename(g_ncid));
	ncAddGlobalAttributes(g_ncid);
	return g_ncid;
}

void ncAddGlobalAttributes(int ncid)
{
	int status; 

	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"ADDA_VERSION",strlen(ADDA_VERSION),ADDA_VERSION))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global attribute 'ADDA_VERSION' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	} 

	time_t now; 
	time(&now); 

	char buf[sizeof "2011-10-08T07:07:09Z"]; // ISO 8601 format
        strftime(buf, sizeof buf, "%FT%TZ", gmtime(&now));

    	char histstrbuf[MAX_LINE];
    	SnprintfErr(ALL_POS,histstrbuf,MAX_LINE,"%s file created by ADDA",buf);

    	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"history",strlen(histstrbuf),histstrbuf))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global attribute 'history' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}
}

void ncFCloseErr(int ncid,ERR_LOC_DECL)
{
	int status;

    	status=nc_close(ncid); 

    	if(status!=NC_NOERR) {
		LogError(ERR_LOC_CALL,"Error closing file '%s', ncerror='%s': ",ncGetFilename(ncid),nc_strerror(status));
	}

    	g_ncmsgbuf[0]=0; // clear fname string buf
}

int ncInitDipFile(int ncid,int *bX,int *bY,int *bZ,int *Nm)
{
	size_t geom_xDim, geom_yDim, geom_zDim;
	int geom_ndims, g_geom_dimids[NC_MAX_VAR_DIMS], geom_natts; // geom_id is global
	char dimensionName[NC_MAX_NAME];
	nc_type geom_type; 

	TIME_TYPE tstart=GET_TIME();

	if (nc_inq_varid(ncid,"geom",&g_geom_id)!=NC_NOERR) {
		LogError(ONE_POS,"No scatterer definition (variable: 'geom') in file %s found",ncGetFilename(ncid));
	}

	if (nc_inq_var(ncid, g_geom_id,0,&geom_type,&geom_ndims,g_geom_dimids,&geom_natts)!=NC_NOERR) {
		LogError(ONE_POS,"Something is wrong with file %s ",ncGetFilename(ncid));
	}

	if (geom_ndims!=3) {
		LogError(ONE_POS,"Wrong dimensions of geom variable, expected 3 got %d file %s ",geom_ndims,ncGetFilename(ncid));
	}

	if (nc_inq_dim(ncid,g_geom_dimids[0],dimensionName,&geom_xDim)!=NC_NOERR) {
		LogError(ONE_POS,"Error reading grid box X size from file %s ",ncGetFilename(ncid));
	}

	if (nc_inq_dim(ncid,g_geom_dimids[1],dimensionName,&geom_yDim)!=NC_NOERR) {
		LogError(ONE_POS,"Error reading grid box Y size from file %s ",ncGetFilename(ncid));
	}

	if (nc_inq_dim(ncid,g_geom_dimids[2],dimensionName,&geom_zDim)!=NC_NOERR) {
		LogError(ONE_POS,"Error reading grid box Z size from file %s ",ncGetFilename(ncid));
	}

	// printf("box size: %zu, %zu, %zu\n", geom_xDim,geom_yDim,geom_zDim);

	if (nc_get_att_int(ncid,g_geom_id,"Nmat",&Nmat)!=NC_NOERR) {
		LogError(ONE_POS,"Error reading Nmat from file %s", ncGetFilename(ncid));
	}

	if (Nmat>MAX_NMAT) {
		LogError(ONE_POS,"Adda was compiled with MAX_NMAT=%d, but we found %d in %s ", MAX_NMAT, Nmat,ncGetFilename(ncid));
	}

	unsigned long long occNdpls;
	if (nc_get_att_ulonglong(ncid,g_geom_id,"N_occupied_dipols",&occNdpls)!=NC_NOERR) {
		LogError(ONE_POS,"Error reading N_occupied_dipols from file %s ",ncGetFilename(ncid));
	} 

	// *rft is set after the call to this function
	*Nm=Nmat;
	*bX=geom_xDim*jagged;
 	*bY=geom_yDim*jagged;
	*bZ=geom_zDim*jagged;

	nvoid_Ndip=occNdpls*jagged*jagged*jagged; 

	// VARPRINTF(nvoid_Ndip); 
	
	Timing_FileIO+=GET_TIME()-tstart;

	return EXIT_SUCCESS;
}

void ncReadDipFile(unsigned char * restrict material_tmp, unsigned short * position_full)
{
	int debugMaterial=0; // debug
	int x0,y0,z0,x,y,z;
	size_t index=0;

#ifndef SPARSE
	// to remove possible overflows
	size_t boxX_l=(size_t)boxX;
#endif 

	TIME_TYPE tstart=GET_TIME(); 
	
        /*
        VARPRINTF(local_z0);
	VARPRINTF(local_z1_coer);
	VARPRINTF(local_Ndip);
	VARPRINTF(boxXY);
        */ 
    
#ifndef SPARSE
	printf("non SPARSE branch of ncReadDipFile\n");
	(void)position_full; // suppress unused var mesg
#else
    	printf("SPARSE branch of ncReadDipFile\n");
#endif // SPARSE

#ifdef PARALLEL
	// ------------------- this loads complete scatterer at once
	unsigned char*geom_tmp=NULL;
    	MALLOC_VECTOR(geom_tmp,uchar,boxX*boxY*boxZ,ALL);

    	if( nc_get_var_ubyte(g_ncid, g_geom_id, geom_tmp) != NC_NOERR ) {
        	LogError(ONE_POS,"No scatterer definition (variable: 'geom') in file %s found",ncGetFilename(g_ncid));
    	}
    	// -----------------

/*
        // --------- read the locally needed z-slabs, or at least in slices
	VARPRINTF(local_z0);
	VARPRINTF(local_z1_coer);
	VARPRINTF(jagged);
	unsigned char*geom_tmp=NULL;
	size_t local_zslab_size = (local_z1_coer/jagged) - (local_z0/jagged);
	
	// MALLOC_VECTOR(geom_tmp,uchar,geom_xDim*geom_yDim*local_zslab_size,ALL);
	MALLOC_VECTOR(geom_tmp,uchar,boxXY*local_zslab_size,ALL);
	
	// define slab
	const size_t start[] = {0, 0, (local_z0/jagged) };
	//const size_t count[] = {geom_xDim, geom_yDim, local_zslab_size}; // edge lengths
	const size_t count[] = {boxX, boxY, local_zslab_size}; // edge lengths
	VARPRINTF(local_zslab_size);
	
	if( nc_get_vara_ubyte (ncid, geom_id, &start, &count, geom_tmp) != NC_NOERR ) {
	LogError(ALL_POS,"No scatterer definition (variable: 'geom') in file %s found",fname);
	}
*/
#else // SEQ
	unsigned char*geom_tmp=NULL;
    	MALLOC_VECTOR(geom_tmp,uchar,boxX*boxY*boxZ,ALL);

    	if( nc_get_var_uchar(g_ncid,g_geom_id, geom_tmp) != NC_NOERR ) {
		LogError(ONE_POS,"No scatterer definition (variable: 'geom') in file %s found",ncGetFilename(g_ncid));
	}
#endif

/*
  	VARPRINTF(boxX);
	VARPRINTF(boxX_l);
	VARPRINTF(boxY);
	VARPRINTF(boxZ);
	VARPRINTF(local_z0);
	VARPRINTF(jagged);
	VARPRINTF(local_Ndip);
*/

    	int mat=0;
    	for(x0=0; x0<boxX; x0++) {
	// printf("x=%u\n", x);
    	for(y0=0; y0<boxY; y0++) {
    	// for(size_t z0=0; z0<geom_zDim; z0++) {
    	for(z0=(local_z0/jagged); z0<(local_z1_coer/jagged); z0++) {

#ifdef PARALLEL
        mat = geom_tmp[ tolinearC3D(x0,y0,z0-((local_z0/jagged)),boxX,boxY,boxZ) ]; 
	// printf("outer(%u,%u,%u)=%u\n",x0,y0,z0,geom_tmp[ tolinearC3D(x0,y0,z0,boxX,boxY,local_zslab_size) ]);
#else
        mat = geom_tmp[ tolinearC3D(x0,y0,z0,boxX,boxY,boxZ) ]; // ok for sphere 
	// printf( "outer(%u,%u,%u)=%u\n", x0,y0,z0,geom_tmp[ tolinearC3D(x0,y0,z0,boxX,boxY,boxZ) ] );
#endif 
	// check for Nmat violation in var geom 
	if (mat>Nmat) {
		LogWarning(EC_WARN,ONE_POS,"Nmat (%d), as given in %s, is not equal to the maximum domain of variable 'geom'.", Nmat,ncGetFilename(g_ncid));
	} 

        if (mat!=0) { // 0 in the nc file is a non occupied dipole
#ifndef SPARSE
		for (z=jagged*z0;z<jagged*(z0+1);z++) {
                if (z>=local_z0 && z<local_z1_coer) {
                    for (x=jagged*x0;x<jagged*(x0+1);x++) {
                    for (y=jagged*y0;y<jagged*(y0+1);y++) {
                        index=(z-local_z0)*boxXY+y*boxX_l+x;
                        // printf("(%d,%d,%d)=[%zu]\n",x,y,z,index);
                            if (material_tmp[index]!=Nmat) {
				// VARPRINTF(material_tmp[index]);
				// VARPRINTF(index);
                                printf("\n");
                                LogError(ONE_POS,"Duplicate dipole was found at line in dipole file %s",ncGetFilename(g_ncid));
                            }

                            material_tmp[index]=(unsigned char)(mat-1);
                    } // for y
                    } // for x
		} // if (local_z0<=z<local_z1_coer)
		} // for z 
#else // SPARSE 
		for (z=0;z<jagged;z++) for (y=0;y<jagged;y++) for (x=0;x<jagged;x++) {
                if ((index >= local_nvoid_d0) && (index < local_nvoid_d1)) {
                    material[index-local_nvoid_d0]=(unsigned char)(mat-1);
                    position_full[3*index]=x0*jagged+x;
                    position_full[3*index+1]=y0*jagged+y;
                    position_full[3*index+2]=z0*jagged+z;
                } 
		index++;
		}
#endif
    	} // if mat!=0, occupied gridbox

	} // for z0
	} // for y0
	} // for x0

	Free_general(geom_tmp);

    	ncFCloseErr(g_ncid,ONE_POS);

    	if(debugMaterial) {
    	for ( z=0; z<boxZ;z++) {
        	printf("z=%u\n", z);
        	if (z>=local_z0 && z<local_z1_coer) {
            	for ( x=0;x<boxX;x++) {
                for ( y=0;y<boxY;y++) {
                    index=(z-local_z0)*boxXY+y*boxX_l+x;
                    printf("%u ", material_tmp[index] );
                    // printf("p(%u,%u,%u) ",position_tmp[index*3],  position_tmp[index*3+1],position_tmp[index*3+2]);
                    // printf( "v %u %u %u \n",position_tmp[index*3],  position_tmp[index*3+1],position_tmp[index*3+2]);
                } // for y
                printf("\n");
            	} // for x
            	printf("\n");
		} // if  (z>=local_z0 && z<local_z1_coer)
	} // for z
	} // if debug

	Timing_FileIO+=GET_TIME()-tstart;
}

void ncSaveGeometry(const char * restrict save_geom_fname)
{
	char fname[MAX_FNAME];

    	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/%s",directory,save_geom_fname); 

	int ncid;
    	ncid=ncFCreateErr(fname,0,ALL_POS);

    	ncAddGlobalAttributes(ncid);

	if (netcdf4_addCFheader){
		ncAddCFAttributes(ncid);
	}

	ncSaveGeometryVar(ncid);

    	ncFCloseErr(ncid,ONE_POS);
}


void ncSaveGeometryVar(int ncid)
{
	int status;

/*
  	printf("ncSaveGeometryVar(  %s\n", ncGetFilename(ncid));
	VARPRINTF(local_z0);
	VARPRINTF(local_z1_coer);
	VARPRINTF(local_Ndip);
	VARPRINTF(boxXY);
	VARPRINTF(boxX);
	VARPRINTF(boxY);
	VARPRINTF(boxZ);
	VARPRINTF(local_nvoid_Ndip);
*/

    	int xdimid, ydimid, zdimid;
    	if ((status=nc_def_dim(ncid,"x",boxX,&xdimid))!=NC_NOERR) {
		LogError(ALL_POS,"Error creating x dimension in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

    	if ((status=nc_def_dim(ncid,"y",boxY,&ydimid))!=NC_NOERR) {
        	LogError(ALL_POS,"Error creating y dimension in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
    	}

    	if ((status=nc_def_dim(ncid, "z",boxZ,&zdimid))!=NC_NOERR) {
        	LogError(ALL_POS,"xdimidError creating y dimension in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
    	} 
    
    	int geom_dimids[3] = {xdimid, ydimid, zdimid};
    	int geom_var_id;

    	if ((status=nc_def_var(ncid,"geom",NC_BYTE,3,geom_dimids,&geom_var_id))!=NC_NOERR) { // we want this, UBYTE is no CF conform
        	LogError(ALL_POS,"Error creating var 'geom' dimension in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

    	if ((status=nc_put_att_int(ncid,geom_var_id,"Nmat",NC_INT,1,&Nmat))!=NC_NOERR) {
       		LogError(ALL_POS,"Error adding attribute 'Nmat' to var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
    	}

    	const int NoccupiedDipoles = nvoid_Ndip;
    	if ((status=nc_put_att_int(ncid, geom_var_id, "N_occupied_dipols", NC_UINT64, 1, &NoccupiedDipoles))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding attribute 'N_occupied_dipols' to var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
    	} 
    
/*
	// netcdf documentation is unclear it seems it is better to use valid_range instead.
	// these make the var CF conform (cfconventions.org/)
	// http://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-f90/Attribute-Conventions.html
	const int valid_min=1; // because fillvalue is set 0, see docs
	const int valid_max=255;
	if ((status=nc_put_att_int(ncid,geom_var_id,"valid_min",NC_INT,1,&valid_min))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding attribute 'valid_min' to var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}
	
	if ((status=nc_put_att_int(ncid,geom_var_id,"valid_max",NC_INT,1,&valid_max))!=NC_NOERR) {
        	LogError(ALL_POS,"Error adding attribute 'valid_max' to var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}
*/

    	int valid_range[2] = {0,255};
    	if ((status=nc_put_att_int(ncid,geom_var_id,"valid_range",NC_INT,2,&valid_range[0]))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding attribute 'valid_range' to var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

    	char _Unsigned[] = "true"; 
    	if ((status=nc_put_att_text(ncid,geom_var_id,"_Unsigned",strlen(_Unsigned), _Unsigned))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global attribute '_Unsigned' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

    	// long_name
    	char long_name[] = "Geometry definition of an electromagnetic scatterer to be used as input to the Amsterdam Discrete Dipole Approximation (ADDA) "\
			"code. This variable defines the scatterer by defining one dipole per voxel along the (unscaled) scatterer Carthesian coordinates "\
			"x,y,z. Different refractive indices (<15 without re-complilation of ADDA) are set per dipole by specifying "\
			"the number out of a refractive index material list given elsewhere. Empty dipoles are set to 0 or better not written (results in"\
			"a NaN entry). NaN can be evaluated to 0 upon reading, as this is the set _FillValue. Data is unsigned.";

    	if ((status=nc_put_att_text(ncid,geom_var_id,"long_name",strlen(long_name),long_name))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'long_name' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
    	} 

    	/*
	// we could use the scale factor, to map to real world size. 
	double scale_factor = 53.2;
	if ((status=nc_put_att_double(ncid,geom_var_id,"scale_factor",NC_DOUBLE,1,&scale_factor))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding attribute 'scale_factor' to var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}
	*/ 
   	/*
	int add_offset = 128;
	if ((status=nc_put_att_int(ncid,geom_var_id,"add_offset",NC_INT,1,&add_offset))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding attribute 'add_offset' to var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}
	*/ 

    	const char units[]="nm"; 
    	// const char units[]="no_unit";  // no_unit is not CF conform
	if ((status=nc_put_att_text(ncid, geom_var_id,"units",strlen(units),units))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'long_name' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	} 

	size_t chunksizes[3]; 
	if ((boxX>=32)&&(boxY>=32)&&(boxZ>=32)) {
		chunksizes[0]=chunksizes[1]=chunksizes[2]=32;
	}
	else {
		chunksizes[0]=boxX; 
		chunksizes[1]=boxY; 
	        chunksizes[2]=boxZ; 
        } 

        // if ((status=nc_def_var_chunking(ncid, geom_var_id, NC_CHUNKED, NULL))!=NC_NOERR) { // is the default, but also compresses same as 32,32,32 but slower 
	if ((status=nc_def_var_chunking(ncid, geom_var_id, NC_CHUNKED, &chunksizes[0]))!=NC_NOERR) { 
       		LogError(ALL_POS,"Error setting chunksize for 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	} 
		
	// default NC_FILL_BYTE = ((signed char)-127), we might need this later 
	//char fillval=-128; // must be same type as var 
	unsigned char fillval=0; 
	int fillmode=NC_FILL; // or, NC_NOFILL which does not work in all cases (garbage in data) 
	if ((status=nc_def_var_fill(ncid, geom_var_id, fillmode, &fillval))!=NC_NOERR) {
		LogError(ALL_POS,"Error setting fillvalue for 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	} 
	
#ifdef PARALLEL
        // parallel netcdf cannot write compressed data 
	// int access=NC_COLLECTIVE; // per var, preferred as lower level can coordinate writes, but caused a deadlock issues once single PC linux (openmpi)
     	int access=NC_INDEPENDENT;
    	if ((status = nc_var_par_access(ncid,geom_var_id,access))!=NC_NOERR) {
		LogError(ALL_POS,"Cannot set NC_COLLECTIVE for var 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
        }
        LogWarning(EC_WARN,ONE_POS,"Cannot use compression for var 'geom' in file '%s' for parallel writes. "\
                            	   "Repack the data using nccopy -d9 -cx/32,y/32,z/32 in.nc out.nc\n",ncGetFilename(ncid));
#else
        int   shuffle = NC_SHUFFLE;
	int   deflate = 1;
        int   deflate_level = 9;
        if ((status=nc_def_var_deflate(ncid, geom_var_id, shuffle, deflate, deflate_level))!=NC_NOERR) {
		LogError(ALL_POS,"Error enabling compression in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
        }
#endif // PARALLEL

    	if ((status=nc_enddef(ncid))!=NC_NOERR) {
        	LogError(ALL_POS,"Error setting deflate_level for 'geom' in file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
    	}

	// write occupied parts of slices (x,y bbox per z-slice only) 
	// fastest running index is in x 

    	unsigned char slicebuf[boxX*boxY]; 
    	unsigned char slicebufreduced[boxX*boxY]; 
	memset(slicebuf, 0, sizeof slicebuf); 
	memset(slicebufreduced, 0, sizeof slicebufreduced); 

/*
  	char slicebuf[boxX*boxY];
  	char slicebufreduced[boxX*boxY];
  	memset(slicebuf, -127, sizeof slicebuf);
  	memset(slicebufreduced, -127, sizeof slicebufreduced);
*/ 

    	size_t sliceminx,sliceminy; // slice bbox 
    	size_t slicemaxx,slicemaxy; 
	sliceminx=sliceminy=SIZE_MAX; 
	slicemaxx=slicemaxy=0; 
	    
    	bool writeslice = false; 
    	size_t j; 
	     // VARPRINTF(local_nvoid_Ndip); 
	for(size_t i=0;i<local_nvoid_Ndip;i++) {
		j=i*3; 
		const size_t xi = position[j]; 
	    	const size_t yi = position[j+1]; 
	    	const size_t zi = position[j+2]; 

		    //printf("pos[%zu]: (%u,%u,%u)=%u\n",i,position[j],position[j+1],position[j+2],material[i]+1); 

		// update the bbox 
		if (xi<sliceminx) sliceminx=xi; 
		if (yi<sliceminy) sliceminy=yi; 

		if (xi>slicemaxx) slicemaxx=xi; 
		if (yi>slicemaxy) slicemaxy=yi; 

		slicebuf[ tolinearC2D(xi,yi,boxX,boxY) ] = (material[i]+1); 

		if (i<(local_nvoid_Ndip-1)) { // look one ahead, if there is one dipole
			const size_t nextzi = position[((i+1)*3)+2]; 
			if (nextzi>zi) { // is the last entry for one z-slice, new z-slice starts afterwards 
				writeslice = true;
		        }    
		} 

		if (i==(local_nvoid_Ndip-1)) { // we are at very last entry in list
			writeslice = true;
		} 

		if (writeslice==true)    
		{
		// re seat the slice into (0,0)
		size_t lx = (slicemaxx-sliceminx)+1; 
		size_t ly = (slicemaxy-sliceminy)+1; 

		size_t ii,jj; 
		for(ii=sliceminx; ii<(sliceminx+lx); ii++) {
		for(jj=sliceminy; jj<(sliceminy+ly); jj++) {
			slicebufreduced[tolinearC2D(ii-sliceminx,jj-sliceminy,lx,ly)] = slicebuf[tolinearC2D(ii,jj,boxX,boxY)]; 
			//printf("%c", slicebufreduced[tolinearC2D(ii-sliceminx,jj-sliceminy,lx,ly)] + '0' );
		} 
		} 
		
		// write the reduced slice 
		const size_t start[] = {sliceminx,sliceminy,zi}; 
		const size_t count[] = {lx,ly,1}; // edge lengths 
			
		if ((status=nc_put_vara_uchar(ncid, geom_var_id,start,count,&slicebufreduced[0]))!=NC_NOERR) {
			LogError(ALL_POS,"error %d, %s writing geom z-slab in file %s found",status,nc_strerror(status),ncGetFilename(ncid));
		} 
		// printf("start[] = {%zu,%zu,%zu}; count[] = {%zu,%zu,%zu}\n",  start[0],start[1],start[2],count[0],count[1],count[2]); 
		
		// empty slice buffers 
		memset(slicebuf, 0, sizeof slicebuf); 
		memset(slicebufreduced, 0, sizeof slicebufreduced); 

		// bbox re init 
		sliceminx=sliceminy=SIZE_MAX; 
		slicemaxx=slicemaxy=0; 

		writeslice=false;
		} 
		
	} // for i
}

void ncSaveComplexField(const complex*field, int geometryinfo)
{
}

void ncAddCFAttributes(int ncid)
{

	// CF attributes, for now hardwired we need a way of communicating them: separate file (simple txt, or nc4) or commandline options.
	int status;

	char Conventions[] = "CF-1.6";
    	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"Conventions",strlen(Conventions),Conventions))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'conventions' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

	char title[] = "Hexagonal ice crystal";
    	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"title",strlen(title),title))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'title' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

	char institution[] = "Institute Affiliation"; 
	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"institution",strlen(institution),institution))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'institution' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

	char source[] = "ADDAorfile e.g."; 
	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"source",strlen(source),source))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'source' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	} 

	char references[] = "Some references."; 
	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"references",strlen(references),references))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'references' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}

	char comment[] = "Testcase of a multiline comment.\nThis is the second line."; 
	if ((status=nc_put_att_text(ncid, NC_GLOBAL,"comment",strlen(comment),comment))!=NC_NOERR) {
		LogError(ALL_POS,"Error adding global CF attribute 'comment' to file '%s', ncerro='%s'",ncGetFilename(ncid),nc_strerror(status));
	}
}

// this returns the file name and puts it into global buffer ncmsgbuf[MAX_FNAME]
const char * ncGetFilename(int ncid)
{
	int status; 
	size_t pathlen;

    	if ((status=nc_inq_path(ncid, &pathlen, NULL))!=NC_NOERR) {
	 	if (status==NC_EBADID) {
			LogError(ALL_POS,"Invalid file handle."); 
			g_ncmsgbuf[0]=0;
		}
		else {
			LogError(ALL_POS,"Error retrieving filename.");
		}
	} 

	if ((pathlen<MAX_FNAME)&&(pathlen>0)) {
		if ((status=nc_inq_path(ncid, &pathlen, g_ncmsgbuf))!=NC_NOERR) {
			LogError(ALL_POS,"Error retrieving filename.");
		}
	}
	else {
		LogError(ALL_POS,"Trying to retrieve a file name longer than MAX_FNAME.");
	} 

	return g_ncmsgbuf;
}

inline size_t  tolinearC2D(size_t i,size_t j,size_t nx,size_t ny )
{
	(void)nx; // suppress unused var mesg
	return  (j + i*ny);
}

inline size_t tolinearC3D (size_t i,size_t j,size_t k,size_t nx,size_t ny,size_t nz)
{
	(void)nx; // suppress unused var mesg 
	return (k + j*nz + i*ny*nz);
}

