#include "FluxResult.h"
#include "../Defintions.h"

namespace Flux
{

CFluxResult::CFluxResult(void)
{
	this->Clear();
}

void CFluxResult::Clear(){
	this->m_flux = 0.0;
	
	// by default, this is a very good measurement
	this->m_fluxQualityFlag = FLUX_QUALITY_GREEN;

	// the errors
	m_fluxError_Wind = 0.0;
	m_fluxError_PlumeHeight = 0.0;

	this->m_windField = CWindField();
	this->m_plumeHeight = CPlumeHeight();
	
	m_numGoodSpectra = 0;

	m_coneAngle = NOT_A_NUMBER;
	m_tilt      = NOT_A_NUMBER;
	m_compass   = NOT_A_NUMBER;
	m_volcano   = -1;
	
	m_startTime    = CDateTime(0,0,0,0,0,0);
	m_stopTime	   = CDateTime(0,0,0,0,0,0);
	m_instrument   = "";
	
	// m_instrumentType	= INSTR_GOTHENBURG;
	m_scanOffset		= 0.0;
	m_completeness		= 0.0;
	m_plumeCentre[0]	= NOT_A_NUMBER;
	m_plumeCentre[1]	= NOT_A_NUMBER;
}

CFluxResult::~CFluxResult(void)
{
}

/** Assignment operator */
CFluxResult &CFluxResult::operator=(const CFluxResult &res){
	this->m_flux	= res.m_flux;
	this->m_fluxQualityFlag = res.m_fluxQualityFlag;

	this->m_fluxError_Wind = res.m_fluxError_Wind;
	this->m_fluxError_PlumeHeight = res.m_fluxError_PlumeHeight;

	this->m_windField = res.m_windField;
	this->m_plumeHeight = res.m_plumeHeight;

	m_numGoodSpectra = res.m_numGoodSpectra;

	m_coneAngle    = res.m_coneAngle;
	m_tilt         = res.m_tilt;
	m_compass      = res.m_compass;
	m_volcano      = res.m_volcano;

	m_startTime    = res.m_startTime;
	m_stopTime	   = res.m_stopTime;

	m_instrument = res.m_instrument;

	// m_instrumentType	= res.m_instrumentType;
	m_scanOffset		= res.m_scanOffset;
	m_completeness		= res.m_completeness;
	m_plumeCentre[0]	= res.m_plumeCentre[0];
	m_plumeCentre[1]	= res.m_plumeCentre[1];


	return *this;
}

}

