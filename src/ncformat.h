/* File: ncformat.h
 * $Date::                            $
 * Descr: NetCDF interface
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
#ifndef __ncformat_h
#define __ncformat_h

#ifdef NETCDF4

#ifdef PARALLEL
    #include <netcdf_par.h>
    // install and download from here: http://www.unidata.ucar.edu/software/netcdf/
    // parallel build: http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
#endif // PARALLEL

#include <netcdf.h>
#include <string.h> //strlen

#include "io.h"  //  ERR_LOC_DECL
#include "vars.h"
#include "memory.h"

// debug util
#define VARPRINTF(v) {printf( #v "= %llu\n", (unsigned long long)(v) );}

extern int g_ncid;
extern char*ncBUF; // used in ncGetFilename;

// 'public:'
int  ncFOpenErr(const char * restrict fname,int mode,ERR_LOC_DECL);
int  ncFCreateErr(const char * restrict fname,int ncmode,ERR_LOC_DECL);
void ncFCloseErr(int ncid, ERR_LOC_DECL);

int  ncInitDipFile(int ncid,int * bX,int * bY,int * bZ,int * NmncReadDipFile);
void ncReadDipFile(unsigned char * restrict material_tmp, unsigned short * position_full);

void ncSaveGeometry(const char * restrict save_geom_fname);
void ncSaveGeometryVar(int ncid);

void ncSaveComplexField(const complex*field, int geometryinfo);

// 'private:'
size_t tolinearC2D(size_t i,size_t j,size_t nx,size_t ny);
size_t tolinearC3D(size_t i,size_t j,size_t k,size_t nx,size_t ny,size_t nz);
void ncAddGlobalAttributes(int ncid);
void ncAddGlobalAttributesCF(int ncid);
void ncAddCFAttributes(int ncid);
const char * ncGetFilename(int lncid);

#endif // NETCDF4

#endif // __ncformat_h
