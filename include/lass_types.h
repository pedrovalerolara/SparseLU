/**
 *
 * 	@file ddss_types.h
 *
 * 	@brief Enumerations.
 *
 * 	LASs is a software package provided by:
 * 	Barcelona Supercomputing Center - Centro Nacional de Supercomputacion
 *
 * 	@author Pedro Valero-Lara pedro.valero@bsc.es
 * 	@date 2017-01-02
 * 	@reviewer 
 * 	@modified 
 *
 **/

#ifndef LASS_TYPES_H
#define LASS_TYPES_H

// Enum. for return.
enum LASS_RETURN {Success=0, NoSuccess=1};

#endif

#ifndef DDSS_TYPES_H
#define DDSS_TYPES_H

// DDSS Enumeration 

// Enum. for row major or column major.
enum DDSS_ORDER {RowMajor=101, ColMajor=102};

// Enum for trans.
enum DDSS_TRANS {NoTrans=111, Trans=112};

// Enum. for uplo.
enum DDSS_UPLO {Upper=121, Lower=122};

// Enum. for diag.
enum DDSS_DIAG {NonUnit=131, Unit=132};

// Enum. for side.
enum DDSS_SIDE {Left=141, Right=142};

#endif
