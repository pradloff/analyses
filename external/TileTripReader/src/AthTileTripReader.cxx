/* 
 * File:   AthTileTripReader.cxx
 * Author: Stephen Cole <stephen.cole@cern.ch>
 * 
 * Created on August 21, 2012, 11:32 AM
 */

#include "AthenaBaseComps/AthAlgorithm.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"
#include "AthenaKernel/errorcheck.h"
#include "EventKernel/INavigable4Momentum.h"
#include "PathResolver/PathResolver.h"
#include "GaudiKernel/MsgStream.h"

#include "TSystem.h"

#include "TileTripReader/AthTileTripReader.h"

using namespace std;
using namespace Root;

AthTileTripReader::AthTileTripReader(const std::string& type, 
        const std::string& name, 
        const IInterface* parent): AthAlgTool(type,name,parent) {
    
    m_tripReader=new Root::TTileTripReader;
    declareInterface<IAthSelectorTool>(this);
    declareInterface<IUserDataCalcTool>(this);
    
    declareProperty("TileTripFile",
            m_tripFile="TileTripReader/CombinedTrips_2011_2012.root",
            "Path to the trip list file");
    
    declareProperty("dR",m_dR=0.4);
        
}


AthTileTripReader::~AthTileTripReader() {
    if(m_tripReader)delete m_tripReader;
}


StatusCode AthTileTripReader::initialize(){
    MsgStream mslog( msgSvc(), name() );
    const char* fname;
    if(m_tripFile.find("/")==0 || m_tripFile.find("$")==0 || 
            m_tripFile.find(".")==0 || m_tripFile.find(":")!=string::npos){
        fname = gSystem->ExpandPathName( m_tripFile.c_str() );
    }
    else{
        m_tripFile = PathResolverFindXMLFile( m_tripFile );
        if ( m_tripFile.empty() ) 
            ATH_MSG_WARNING("Could NOT resolve file name "<<m_tripFile);
        fname=m_tripFile.data();
    }
    
    if(!m_tripReader->setTripFile(fname)){
        ATH_MSG_ERROR("ERROR! No trip file loaded!");
        return StatusCode::FAILURE;
    }
    MSG::Level cLevel=mslog.currentLevel();
    if(cLevel<=MSG::DEBUG)
        m_tripReader->setVerbosity(1);
    else m_tripReader->setVerbosity(0);
    m_tripReader->initialize();
    return StatusCode::SUCCESS;
}

StatusCode AthTileTripReader::finalize(){
    if(!m_tripReader->finalize()){
        ATH_MSG_ERROR ( "ERROR! Something went wrong at finalize!" );
        return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
}

const TAccept& AthTileTripReader::accept(const INavigable4Momentum* part){
    const EventInfo* evt=0;
    StatusCode sc=evtStore()->retrieve(evt,"ByteStreamEventInfo");
    if(sc.isFailure()){
        ATH_MSG_ERROR ( "ERROR! Could not retrieve event info from store gate!");
        return m_acceptDummy;
    }
    const EventID& id=*evt->event_ID();
    return m_tripReader->accept(id.run_number(),
            id.lumi_block(),
            id.event_number(),
            part->eta(),
            part->phi());
}

const TResult& AthTileTripReader::calculate(const INavigable4Momentum* part){
    const EventInfo* evt=0;
    StatusCode sc=evtStore()->retrieve(evt,"ByteStreamEventInfo");
    if(sc.isFailure()){
        ATH_MSG_ERROR ( "ERROR! Could not retrieve event info from store gate!");
        return m_resultDummy;
    }
    const EventID& id=*evt->event_ID();    
    return m_tripReader->calculate(id.run_number(),
            id.lumi_block(),
            part->eta(),
            part->phi(),
            m_dR);
}

StatusCode AthTileTripReader::calculateElementUserData(const IAthenaBarCode* abc){
    // Test if we got a jet
    const INavigable4Momentum* inav = dynamic_cast<const INavigable4Momentum*>(abc);
    if ( !inav ){
        ATH_MSG_ERROR ( "Could NOT cast to a const INavigable4Momentum!!!" );
        return StatusCode::FAILURE;
    }

    // get the result of the selection
    const Root::TAccept& myAccept = this->accept(inav);
    bool isAccepted = bool(myAccept);
    ATH_MSG_VERBOSE( "This one is accepted: " << isAccepted );

    // Store the result into user data (the name of the user data should be configurable)
    ATH_CHECK ( userStore()->record(*inav, "passTileTripCleaning", isAccepted) );

    // get the result of the fraction
    const Root::TResult& myResult = this->calculate(inav);
    double myFraction = myResult.getResult("TripAreaFrac");
    ATH_MSG_VERBOSE( "This one has a fractional area with TileTrip of: " << myFraction );

    // Store the result into user data (the name of the user data should be configurable)
    ATH_CHECK ( userStore()->record(*inav, "fractionAreaInTileTrip", myFraction) );

    return StatusCode::SUCCESS;
}
