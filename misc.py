from common.functions import event_function
import ROOT

class count_primary_vertices(event_function):
	def __init__(self,*args,**kwargs):
		event_function.__init__(self,*args,**kwargs)

		self.required_branches+= [
			'vxp_nTracks',
			]
		self.create_branches.update(dict((branch_name,branch_type) for branch_name,branch_type in [
			('nPV_2trks','int'),
			('nPV_3trks','int'),
			('nPV_4trks','int'),
			]))

	def __call__(self,event):
		for i in [2,3,4]:
			event.__dict__['nPV_{0}trks'.format(i)] = sum(1 for nTracks in event.vxp_nTracks if nTracks>=i)

def vector_attributes(collection,attributes):
	d = {}
	for attribute,type_ in attributes.items(): 
		d[attribute] = ROOT.std.vector(type_)()
	for key,value in sorted(collection.items()):
		for attribute in attributes.keys():
			d[attribute].push_back(getattr(value,attribute))
	return d


def list_attributes(collection,attributes,collection_name):
	d = {}
	for attribute in attributes: 
		d[collection_name+attribute] = []
	for key,value in sorted(collection.items()):
		for attribute in attributes:
			d[collection_name+attribute].append(getattr(value,attribute))
	return d

"""
bool hotTileVeto(){

  bool hotTileVeto = false;

  // this veto is ONLY applied to data and NOT to MC

  if (isData || isEmbedding){
    // adding procedure for cleaning
    // Hot Tile calorimeter in period B1 and B2
    if ( RunNumber == 202660 ||
         RunNumber == 202668 ||
         RunNumber == 202712 ||
         RunNumber == 202740 ||
         RunNumber == 202965 ||
         RunNumber == 202987 ||
         RunNumber == 202991 ||
         RunNumber == 203027 ||
         RunNumber == 203169) {
      for(std::vector<TLorentzVector>::size_type i=0; i<theJets.size(); i++) {
        TLorentzVector jet = theJets[i];
        if(jet.Pt()/GeV>20) {
          Float_t j_fmax= jet_AntiKt4TopoEM_fracSamplingMax;
          Int_t   j_smax= jet_AntiKt4TopoEM_SamplingMax;
          Float_t j_eta = jet_AntiKt4TopoEM_eta;
          Float_t j_phi = jet_AntiKt4TopoEM_phi;

          bool _etaphi28 = false;
          if(j_eta>-0.2 && j_eta<-0.1 && j_phi>2.65 && j_phi< 2.75) _etaphi28=true;
          if(j_fmax>0.6 && j_smax==13 && _etaphi28) allJetClean = false;
      }
    }
  }
  return hotTileVeto;
}

"""
