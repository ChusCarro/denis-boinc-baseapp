/*======================================================================
                                                                        
                      --- DENIS Project ---                             
                        -----------------                               
                                                                        
      Distributed computing                                             
      Electrophysiologycal Models                                       
      Networking colaboration                                           
      In Silico research                                                
      Sharing Knowledge                                                 
                                                                        
------------------------------------------------------------------------
                                                                        
 -- C Models Header created for the DENIS Project --                    
                                                                        
 This file has been automatically created using the CellML API.         
                                                                        
------------------------------------------------------------------------
                                                                        
 DENIS-CellML2C Copyright 2015 J. Carro; J. Castro                      
                                                                        
 Licensed under the Apache License, Version 2.0 (the "License");      
 you may not use this file except in compliance with the License.       
 You may obtain a copy of the License at                                
                                                                        
 http://www.apache.org/licenses/LICENSE-2.0                             
                                                                        
 Unless required by applicable law or agreed to in writing, software    
 distributed under the License is distributed on an "AS IS" BASIS,    
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or        
 implied. See the License for the specific language governing           
 permissions and limitations under the License.                         
                                                                        
------------------------------------------------------------------------
                                                                        
                                                  San Jorge University  
                                                 School of Engineering  
                                                                        
                                                     http://eps.usj.es  
                                                                        
======================================================================*/

#include "ten_tusscher_model_2006_IK1Ko_epi_units.h"
#include "tentusscher_noble_noble_panfilov_2004_a.h"
#include "grandi_pasqualini_bers_2010_ss.h"
#include "ten_tusscher_model_2006_IK1Ko_endo_units.h"
#include "tentusscher_noble_noble_panfilov_2004_c.h"
#include "ten_tusscher_model_2006_IK1Ko_M_units.h"
#include "CRLP2011_EPI.h"
#include "tentusscher_noble_noble_panfilov_2004_b.h"
#include "Ohara_Rudy_2011.h"

#define ten_tusscher_model_2006_IK1Ko_epi_units_ID 0
#define tentusscher_noble_noble_panfilov_2004_a_ID 1
#define grandi_pasqualini_bers_2010_ss_ID 2
#define ten_tusscher_model_2006_IK1Ko_endo_units_ID 3
#define tentusscher_noble_noble_panfilov_2004_c_ID 4
#define ten_tusscher_model_2006_IK1Ko_M_units_ID 5
#define CRLP2011_EPI_ID 6
#define tentusscher_noble_noble_panfilov_2004_b_ID 7
#define Ohara_Rudy_2011_ID 8

int getModelId(const char* modelName){
	if(strcmp(modelName, "ten_tusscher_model_2006_IK1Ko_epi_units")==0)
		return ten_tusscher_model_2006_IK1Ko_epi_units_ID;
	if(strcmp(modelName, "tentusscher_noble_noble_panfilov_2004_a")==0)
		return tentusscher_noble_noble_panfilov_2004_a_ID;
	if(strcmp(modelName, "grandi_pasqualini_bers_2010_ss")==0)
		return grandi_pasqualini_bers_2010_ss_ID;
	if(strcmp(modelName, "ten_tusscher_model_2006_IK1Ko_endo_units")==0)
		return ten_tusscher_model_2006_IK1Ko_endo_units_ID;
	if(strcmp(modelName, "tentusscher_noble_noble_panfilov_2004_c")==0)
		return tentusscher_noble_noble_panfilov_2004_c_ID;
	if(strcmp(modelName, "ten_tusscher_model_2006_IK1Ko_M_units")==0)
		return ten_tusscher_model_2006_IK1Ko_M_units_ID;
	if(strcmp(modelName, "CRLP2011_EPI")==0)
		return CRLP2011_EPI_ID;
	if(strcmp(modelName, "tentusscher_noble_noble_panfilov_2004_b")==0)
		return tentusscher_noble_noble_panfilov_2004_b_ID;
	if(strcmp(modelName, "Ohara_Rudy_2011")==0)
		return Ohara_Rudy_2011_ID;
	return -1;
}

int getNameId(int modelId, const char* variable, const char* component, const char** names, int namesLength){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::getNameId(variable, component, names, namesLength);
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::getNameId(variable, component, names, namesLength);
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::getNameId(variable, component, names, namesLength);
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::getNameId(variable, component, names, namesLength);
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::getNameId(variable, component, names, namesLength);
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::getNameId(variable, component, names, namesLength);
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::getNameId(variable, component, names, namesLength);
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::getNameId(variable, component, names, namesLength);
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::getNameId(variable, component, names, namesLength);
	default:
		fprintf(stderr, "Model not implemented\n");
		return -1;
	}
}

double* getNewRatesArray(int modelId){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::getNewRatesArray();
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::getNewRatesArray();
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::getNewRatesArray();
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::getNewRatesArray();
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::getNewRatesArray();
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::getNewRatesArray();
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::getNewRatesArray();
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::getNewRatesArray();
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::getNewRatesArray();
	default:
		fprintf(stderr, "Model not implemented\n");
		return NULL;
	}
}

double* getNewStatesArray(int modelId){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::getNewStatesArray();
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::getNewStatesArray();
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::getNewStatesArray();
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::getNewStatesArray();
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::getNewStatesArray();
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::getNewStatesArray();
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::getNewStatesArray();
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::getNewStatesArray();
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::getNewStatesArray();
	default:
		fprintf(stderr, "Model not implemented\n");
		return NULL;
	}
}

double* getNewAlgebraicArray(int modelId){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::getNewAlgebraicArray();
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::getNewAlgebraicArray();
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::getNewAlgebraicArray();
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::getNewAlgebraicArray();
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::getNewAlgebraicArray();
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::getNewAlgebraicArray();
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::getNewAlgebraicArray();
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::getNewAlgebraicArray();
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::getNewAlgebraicArray();
	default:
		fprintf(stderr, "Model not implemented\n");
		return NULL;
	}
}

double* getNewConstantsArray(int modelId){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::getNewConstantsArray();
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::getNewConstantsArray();
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::getNewConstantsArray();
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::getNewConstantsArray();
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::getNewConstantsArray();
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::getNewConstantsArray();
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::getNewConstantsArray();
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::getNewConstantsArray();
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::getNewConstantsArray();
	default:
		fprintf(stderr, "Model not implemented\n");
		return NULL;
	}
}

void getNewArrays(int modelId, double **constants, double **rates, double **states, double **algebraic){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::getNewArrays(constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::getNewArrays(constants, rates, states, algebraic);
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::getNewArrays(constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::getNewArrays(constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::getNewArrays(constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::getNewArrays(constants, rates, states, algebraic);
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::getNewArrays(constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::getNewArrays(constants, rates, states, algebraic);
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::getNewArrays(constants, rates, states, algebraic);
	default:
		fprintf(stderr, "Model not implemented\n");
	}
}

void getVectorsLength(int modelId, int* algebraicLength, int* statesLength, int* constantsLength){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::getVectorsLength(algebraicLength, statesLength, constantsLength);
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::getVectorsLength(algebraicLength, statesLength, constantsLength);
	default:
		fprintf(stderr, "Model not implemented\n");
	}
}

void names(int modelId, const char* &VoI, const char** constants, const char** rates, const char** states, const char** algebraic){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::names(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::names(VoI, constants, rates, states, algebraic);
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::names(VoI, constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::names(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::names(VoI, constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::names(VoI, constants, rates, states, algebraic);
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::names(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::names(VoI, constants, rates, states, algebraic);
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::names(VoI, constants, rates, states, algebraic);
	default:
		fprintf(stderr, "Model not implemented\n");
	}
}

void initConsts(int modelId, double* constants, double* states){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::initConsts(constants, states);
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::initConsts(constants, states);
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::initConsts(constants, states);
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::initConsts(constants, states);
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::initConsts(constants, states);
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::initConsts(constants, states);
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::initConsts(constants, states);
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::initConsts(constants, states);
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::initConsts(constants, states);
	default:
		fprintf(stderr, "Model not implemented\n");
	}
}

void computeRates(int modelId, double VoI, double* constants, double* rates, double* states, double* algebraic){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::computeRates(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::computeRates(VoI, constants, rates, states, algebraic);
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::computeRates(VoI, constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::computeRates(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::computeRates(VoI, constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::computeRates(VoI, constants, rates, states, algebraic);
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::computeRates(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::computeRates(VoI, constants, rates, states, algebraic);
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::computeRates(VoI, constants, rates, states, algebraic);
	default:
		fprintf(stderr, "Model not implemented\n");
	}
}

void computeVariables(int modelId, double VoI, double* constants, double* rates, double* states, double* algebraic){
	switch(modelId){
	case ten_tusscher_model_2006_IK1Ko_epi_units_ID:
		return ten_tusscher_model_2006_IK1Ko_epi_units::computeVariables(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_a_ID:
		return tentusscher_noble_noble_panfilov_2004_a::computeVariables(VoI, constants, rates, states, algebraic);
	case grandi_pasqualini_bers_2010_ss_ID:
		return grandi_pasqualini_bers_2010_ss::computeVariables(VoI, constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_endo_units_ID:
		return ten_tusscher_model_2006_IK1Ko_endo_units::computeVariables(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_c_ID:
		return tentusscher_noble_noble_panfilov_2004_c::computeVariables(VoI, constants, rates, states, algebraic);
	case ten_tusscher_model_2006_IK1Ko_M_units_ID:
		return ten_tusscher_model_2006_IK1Ko_M_units::computeVariables(VoI, constants, rates, states, algebraic);
	case CRLP2011_EPI_ID:
		return CRLP2011_EPI::computeVariables(VoI, constants, rates, states, algebraic);
	case tentusscher_noble_noble_panfilov_2004_b_ID:
		return tentusscher_noble_noble_panfilov_2004_b::computeVariables(VoI, constants, rates, states, algebraic);
	case Ohara_Rudy_2011_ID:
		return Ohara_Rudy_2011::computeVariables(VoI, constants, rates, states, algebraic);
	default:
		fprintf(stderr, "Model not implemented\n");
	}
}

