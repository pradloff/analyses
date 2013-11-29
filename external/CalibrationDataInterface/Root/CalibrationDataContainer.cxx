///////////////////////////////////////////////////////////////////
// CalibrationDataContainer.cxx, (c) ATLAS Detector software
///////////////////////////////////////////////////////////////////

#include "CalibrationDataInterface/CalibrationDataContainer.h"

#include <iostream>
#include <cassert>
#include <limits>

#include "TH1.h"
#include "TAxis.h"
#include "TF1.h"
#include "TMatrixT.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TString.h"
#include "TObjString.h"

using Analysis::UncertaintyResult;
using Analysis::CalibrationDataContainer;
using Analysis::CalibrationDataHistogramContainer;
using Analysis::CalibrationDataMappedHistogramContainer;
using Analysis::CalibrationDataFunctionContainer;

// Things that are best hidden from the outside world...
namespace {

  // The value below is added to (subtracted from) the lower (upper) bound of validity
  // when a variable is restricted to be within its validity range.
  const double rangeEpsilon = 1.e-5;

  // array size for boundary specifications
  const int maxParameters = 10;

}




// ---------------------- CalibrationDataContainer base class ----------------------------

CalibrationDataContainer::CalibrationDataContainer(const char* name) :
  TMap(), m_objResult(0), m_objSystematics(0), m_restrict(false)
{
  SetName(name);
}

CalibrationDataContainer::~CalibrationDataContainer()
{
}

// short-hand for the total systematic uncertainty retrieval

CalibrationDataContainer::CalibrationStatus
CalibrationDataContainer::getSystUncertainty(const CalibrationDataVariables& x,
					     UncertaintyResult& result, TObject* obj) const
{
  // cache the pointer to the "systematics" object (to avoid string comparisons)
  if (!obj) {
    if (! m_objSystematics) {
      // std::cout << "retrieving total systematics pointer" << std::endl;
      m_objSystematics = GetValue("systematics");
    }
    obj = m_objSystematics;
  }
  return getUncertainty("systematics", x, result, obj);
}

// retrieve the list of uncertainties for this calibration

std::vector<std::string>
CalibrationDataContainer::listUncertainties() const
{
  std::vector<std::string> uncertainties;
  TIter it(GetTable());
  while (TPair* pair = (TPair*) it()) {
    std::string spec(pair->Key()->GetName());
    uncertainties.push_back(spec);
  }
  return uncertainties;
}

// retrieve all uncertainties for this calibration

CalibrationDataContainer::CalibrationStatus
CalibrationDataContainer::getUncertainties(const CalibrationDataVariables& x,
					   std::map<std::string, UncertaintyResult>& all) const
{
  CalibrationStatus mycode(kSuccess);
  UncertaintyResult result;

  // first treat the "result" entry separately
  double single_result;
  CalibrationStatus code = getResult(x, single_result);
  if (code == kError) {
    std::cerr << "in CalibrationDataContainer::getUncertainties(): error retrieving result!" << std::endl;
    return code;
  }
  else if (code != kSuccess) mycode = code;
  result.first = single_result;
  result.second = 0;
  all[std::string("result")] = result;

  // similar for the "statistics" entry
  code = getStatUncertainty(x, single_result);
  if (code == kError) {
    std::cerr << "in CalibrationDataContainer::getUncertainties(): error retrieving stat. uncertainty!" << std::endl;
    return code;
  }
  else if (code != kSuccess) mycode = code;
  result.first  =  single_result;
  result.second = -single_result;
  all[std::string("statistics")] = result;

  // then cycle through the other (systematic) uncertainties
  TIter it(GetTable());
  while (TPair* pair = (TPair*) it()) {
    std::string spec(pair->Key()->GetName());
    // ignore these specific entries
    if (spec == "comment" || spec == "result" || spec == "statistics") continue;
    code = getUncertainty(spec, x, result, pair->Value());
    // we should never be finding any errors
    if (code == kError) {
      std::cerr << "in CalibrationDataContainer::getUncertainties(): error retrieving named uncertainty "
		<< spec << "!" << std::endl;
      return code;
    }
    // this assumes that non-success codes are likely to be correlated between uncertainty sources
    else if (code != kSuccess) mycode = code;
    all[spec] = result;
  }
  return mycode;
}

// retrieve the comments for this calibration (if any)

std::string
CalibrationDataContainer::getComment() const
{
  TObject* obj = GetValue("comment");
  if (! obj) return std::string("");
  TObjString* s = dynamic_cast<TObjString*>(obj);
  if (! s ) return std::string("");
  return std::string(s->GetName());
}

// insert the given object at the position indicated by the given 'uncertainty' index

void
CalibrationDataContainer::setUncertainty(const std::string& unc, TObject* obj)
{
  if (TPair* p = (TPair*) FindObject(unc.c_str())) DeleteEntry(p->Key());
  Add(new TObjString(unc.c_str()), obj);
}

// specialization of the above: insert the calibration result

void
CalibrationDataContainer::setResult(TObject* obj)
{
  setUncertainty(std::string("result"), obj);
}

// special case: TString itself doesn't inherit from TObject

void
CalibrationDataContainer::setComment(const std::string& text)
{
  if (TPair* p = (TPair*) FindObject("comment")) DeleteEntry(p->Key());
  Add(new TObjString("comment"), new TObjString(text.c_str()));
}

int
CalibrationDataContainer::typeFromString(const std::string& key) const
{
  if      (key == "eta") return kEta;
  else if (key == "abseta") return kAbsEta;
  else if (key == "pt") return kPt;
  else if (key == "tagweight") return kTagWeight;
  // return value for unknown keywords
  else return -1;
}


// Determine which variables are to be used, and insert them in a separate array.
// Identical functionality as the preceding function but specialised to the "results" object
// The return value is used to indicate whether any input co-ordinate was out of bounds

bool
CalibrationDataContainer::computeVariables(const CalibrationDataVariables& x) const
{
  // ensure that the variable types have been computed properly
  if (m_variables.size() == 0) computeVariableTypes();

  // also keep track of whether the variables are within bounds
  bool OK(true);

  // for (unsigned int var = 0; var < m_variablesResult->size(); ++var) {
  //   switch ((*m_variablesResult)[var]) {
  for (unsigned int var = 0; var < m_variables.size(); ++var) {
    switch (m_variables[var]) {
    case kPt:
      // assume that the input values are given in MeV but the performance calibration in GeV!
      m_vars[var] = x.jetPt * 0.001;
      break;
    case kEta:
      m_vars[var] = x.jetEta;
      break;
    case kAbsEta:
      m_vars[var] = x.jetEta;
      if (m_vars[var] < 0) m_vars[var] *= -1.0;
      break;
    case kTagWeight:
      m_vars[var] = x.jetTagWeight;
    }
    if (m_vars[var] < getLowerBound(m_variables[var])) {
      OK = false;
      if (m_restrict) m_vars[var] = getLowerBound(m_variables[var]) + rangeEpsilon;
    } else if (m_vars[var] >= getUpperBound(m_variables[var])) {
      OK = false;
      if (m_restrict) m_vars[var] = getUpperBound(m_variables[var]) - rangeEpsilon;
    }
  }

  return OK;
}

double
CalibrationDataContainer::getLowerBound(unsigned int vartype) const
{
  double minDefault = (vartype == kAbsEta || vartype == kPt) ? 0 : -std::numeric_limits<double>::max();
  return (vartype < m_lowerBounds.size()) ? m_lowerBounds[vartype] : minDefault;
}

double
CalibrationDataContainer::getUpperBound(unsigned int vartype) const
{
  return (vartype <  m_lowerBounds.size()) ? m_upperBounds[vartype] : std::numeric_limits<double>::max();
}

std::vector<std::pair<double, double> >
CalibrationDataContainer::getBounds() const
{
  // ensure that the variable types have been computed properly
  // if (m_variables.size() == 0) computeVariableTypes();
  if (m_variables.size() == 0) computeVariableTypes();

  std::vector<std::pair<double, double> > bounds;
  for (unsigned int t = 0; t < m_lowerBounds.size() && t <= kAbsEta; ++t) {
    bounds.push_back(std::make_pair<double, double>(m_lowerBounds[t], m_upperBounds[t]));
  }
  return bounds;
}

// ---------------------- CalibrationDataHistogramContainer class ----------------------

CalibrationDataHistogramContainer::CalibrationDataHistogramContainer(const char* name) :
  CalibrationDataContainer(name), m_interpolate(false)
{
  // Reset the validity bounds to reflect 'no bounds'.
  // They will be re-determined upon the first computation.

  m_lowerBounds.clear();
  m_lowerBounds.resize(maxParameters, -std::numeric_limits<double>::max());
  m_lowerBounds[kPt] = m_lowerBounds[kAbsEta] = 0;
  m_upperBounds.clear();
  m_upperBounds.resize(maxParameters, std::numeric_limits<double>::max());

  // But by default, switch on the range checking
  restrictToRange(true);
}

CalibrationDataHistogramContainer::~CalibrationDataHistogramContainer()
{
}


void
CalibrationDataHistogramContainer::computeVariableTypes() const
{
  // cache pointer to central values histogram
  if (! m_objResult) m_objResult = GetValue("result");

  // histograms need a special treatment, as the coordinate titles are not actually stored
  // with the title itself, but instead moved off to the axis titles...
  const TH1* hobj = dynamic_cast<const TH1*>(m_objResult);

  // no protection against null pointers here -- should not be necessary?
  int dims = hobj->GetDimension();
  for (int dim = 0; dim < dims; ++dim) {
    TAxis* axis;
    switch (dim) {
    case 0: axis = hobj->GetXaxis(); break;
    case 1: axis = hobj->GetYaxis(); break;
    default: axis = hobj->GetZaxis();
    }
    int vartype = typeFromString(axis->GetTitle());
    if (vartype < 0) {
      // Only flag the issue but otherwise take no action (assume non-argument use of a semicolon)
      std::cerr << "in CalibrationDataHistogramContainer::computeVariableTypes(): cannot construct variable type from name "
		<< axis->GetTitle() << std::endl;
    } else {
      m_variables.push_back((unsigned int) vartype);
    }
  }

  // After doing this, we should always have a non-null vector!
  assert(m_variables.size() > 0);

  // Also compute the validity bounds for this calibration object
  const_cast<CalibrationDataHistogramContainer*>(this)->checkBounds();
}
// result retrieval

CalibrationDataContainer::CalibrationStatus
CalibrationDataHistogramContainer::getResult(const CalibrationDataVariables& x,
					     double& result, TObject* obj) const
{
  if (!obj) {
    if (! m_objResult) {
      // std::cout << "retrieving central value pointer" << std::endl;
      m_objResult = GetValue("result");
    }
    obj = m_objResult;
  }
  TH1* hist = dynamic_cast<TH1*>(obj);
  if (! hist) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);
  // find the relevant "global" bin.
  // Note the limitation: at most three dimensions are supported.
  // TH1::FindFixBin() will ignore the variables not needed.
  // Note: FindFixBin() is only available in "recent" ROOT versions (FindBin() is appropriate for older versions)
  // (otherwise we need to rely on the ResetBit(TH1::kCanRebin) method having been used)
  if (m_interpolate) {
    // std::cout << "retrieving interpolated result" << std::endl;
    result = getInterpolatedResult(hist);
  } else {
    // std::cout << "retrieving binned result" << std::endl;
    int bin = hist->FindFixBin(m_vars[0], m_vars[1], m_vars[2]);
    result = hist->GetBinContent(bin);
  }

  return inRange ? kSuccess : kRange;
}

// statistical uncertainty retrieval (special since it is stored with the result itself)

CalibrationDataContainer::CalibrationStatus
CalibrationDataHistogramContainer::getStatUncertainty(const CalibrationDataVariables& x,
						      double& result) const
{
  if (! m_objResult) {
    // std::cout << "retrieving central value pointer" << std::endl;
    m_objResult = GetValue("result");
  }
  TH1* hist = dynamic_cast<TH1*>(m_objResult);
  if (! hist) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);
  // find the relevant "global" bin.
  // Note the limitation: at most three dimensions are supported.
  // TH1::FindFixBin() will ignore the variables not needed.
  // Note: FindFixBin() is only available in "recent" ROOT versions (FindBin() is appropriate for older versions)
  // (otherwise we need to rely on the ResetBit(TH1::kCanRebin) method having been used)
  if (m_interpolate) {
    // interpolating the uncertainties doesn't seem very sensible..
    result = 0;
  } else {
    int bin = hist->FindFixBin(m_vars[0], m_vars[1], m_vars[2]);
    result = hist->GetBinError(bin);
  }

  return inRange ? kSuccess : kRange;
}

// general uncertainty retrieval

CalibrationDataContainer::CalibrationStatus
CalibrationDataHistogramContainer::getUncertainty(const std::string& unc,
						  const CalibrationDataVariables& x,
						  UncertaintyResult& result, TObject* obj) const
{
  // treat statistical uncertainties separately (they are stored with the actual result)
  if (unc == "statistics") {
    double res;
    CalibrationStatus code = getStatUncertainty(x, res);
    if (code == kError) return code;
    result.first =   res;
    result.second = -res;
    return code;
  }

  if (!obj) obj = GetValue(unc.c_str());
  TH1* hist = dynamic_cast<TH1*>(obj);
  if (! hist) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);

  if (m_interpolate) {
    result.first = getInterpolatedResult(hist);
    // symmetrise the uncertainty (as there is no code to interpolate the bin errors)
    result.second = -result.first;
  } else {
    // TH1::FindFixBin() will ignore the variables not needed.
    // Note: FindFixBin() is only available in "recent" ROOT versions (FindBin() is appropriate for older versions)
    // (otherwise we need to rely on the ResetBit(TH1::kCanRebin) method having been used)
    int bin = hist->FindFixBin(m_vars[0], m_vars[1], m_vars[2]);
    // the "first" and "second" entries are filled with the
    // "positive" and "negative" uncertainties, respectively.
    result.first = hist->GetBinContent(bin);
    result.second = hist->GetBinError(bin);
  }

  return inRange ? kSuccess : kRange;
}

// check the bounds of validity for this calibration object

void
CalibrationDataHistogramContainer::checkBounds()
{
  const TH1* hist = dynamic_cast<const TH1*>(m_objResult);
  if (!hist) {
    std::cerr << "in CalibrationDataHistogramContainer::checkBounds(): object type does not derive from TH1" << std::endl;
    return;
  } else if (hist->GetDimension() != int(m_variables.size())) {
    std::cerr << "in CalibrationDataHistogramContainer::checkBounds(): given number of variable types ("
	      << m_variables.size() << ") doesn't match histogram dimension ("
	      << hist->GetDimension() << ")!" << std::endl;
    return;
  }
  for (unsigned int t = 0; int(t) < hist->GetDimension(); ++t) {
    TAxis* axis;
    switch (t) {
    case 0: axis = hist->GetXaxis(); break;
    case 1: axis = hist->GetYaxis(); break;
    default: axis = hist->GetZaxis();
    }
    
    // for (unsigned int t = 0; t < m_variables.size(); ++t) {
    //   if (m_variables[t] > m_upperBounds.size()) {
    // 	std::cerr << "in CalibrationDataHistogramContainer::checkBounds(): variable " << t << "type ("
    // 		  << m_variables[t] << "exceeds maximum type number (" << m_upperBounds.size() << ")!"
    // 		  << std::endl;
    // 	return;
    //   }
    // }
    double amax = axis->GetXmax(), amin = axis->GetXmin();
    if (amax < m_upperBounds[m_variables[t]]) m_upperBounds[m_variables[t]] = amax;
    if (amin > m_lowerBounds[m_variables[t]]) m_lowerBounds[m_variables[t]] = amin;
  }
}

/** indicate whether the given uncertainty is correlated from bin to bin or not
    (note that this function is to be used only for _systematic_ uncertainties) */
bool
CalibrationDataHistogramContainer::isBinCorrelated(const std::string& unc) const
{
  return (m_uncorrelatedSyst.FindObject(unc.c_str()) == 0);
}

/** indicate that the given uncertainty is to be treated uncorrelated from bin to bin
    (note that the default is for all systematic uncertainties to be treated as correlated) */
void
CalibrationDataHistogramContainer::setUncorrelated(const std::string& unc)
{
  m_uncorrelatedSyst.Add(new TObjString(unc.c_str()));
}

/** Indicate whether results are to be interpolated between bins or not
    (this feature is thought to be useful mostly for MC efficiencies) */
void
CalibrationDataHistogramContainer::setInterpolated(bool doInterpolate)
{
  m_interpolate = doInterpolate;
}

/** Indicate whether histogram interpolation is used or not */
bool
CalibrationDataHistogramContainer::isInterpolated() const
{
  return m_interpolate;
}

/** Retrieve interpolated result (utility function) */
double
CalibrationDataHistogramContainer::getInterpolatedResult(TH1* hist) const
{
  switch (hist->GetDimension()) {
  case 3:
    return hist->Interpolate(m_vars[0], m_vars[1], m_vars[2]);
  case 2:
    return hist->Interpolate(m_vars[0], m_vars[1]);
  case 1:
  default:
    return hist->Interpolate(m_vars[0]);
  }
}


/** test whether this calibration object is one for "continuous" calibration
    (this has some subtle consequences for the treatment of bin-to-bin correlations).
    The return value will be -1 in case this is not a "continuous" calibration object,
    and the axis number (0 for X, 1 for Y, 2 for Z) otherwise. */
int
CalibrationDataHistogramContainer::getTagWeightAxis() const
{
  for (unsigned int type = 0; type < m_variables.size(); ++type)
    if (m_variables[type] == kTagWeight) return int(type);
  return -1;
}


// ---------------------- CalibrationDataMappedHistogramContainer class ----------------

CalibrationDataMappedHistogramContainer::CalibrationDataMappedHistogramContainer(const char* name) :
  CalibrationDataHistogramContainer(name), m_lastBin(0)
{
}

CalibrationDataMappedHistogramContainer::~CalibrationDataMappedHistogramContainer()
{
}

void
CalibrationDataMappedHistogramContainer::computeVariableTypes() const
{
  // cache pointer to central values histogram
  if (! m_objResult) m_objResult = GetValue("result");

  // histograms need a special treatment, as the coordinate titles are not actually stored
  // with the title itself, but instead moved off to the axis titles...
  const TH1* hobj = dynamic_cast<const TH1*>(m_objResult);

  // no protection against null pointers here -- should not be necessary?
  int dims = hobj->GetDimension();
  for (int dim = 0; dim < dims; ++dim) {
    TAxis* axis;
    switch (dim) {
    case 0: axis = hobj->GetXaxis(); break;
    case 1: axis = hobj->GetYaxis(); break;
    default: axis = hobj->GetZaxis();
    }
    std::string var(axis->GetTitle());
    if (var == "mapped") {
      // Special case: mapped variables, so make sure to specify the original variables (not the mapped ones).
      // Note that the code here assumes that the mapping is identical for all objects..
      for (unsigned int m = 0; m < m_mapped.size(); ++m) {
	int vartype = typeFromString(m_mapped[m]);
	// this check should never fail; therefore, bail out if this does happen
	assert (! (vartype < 0));
	m_variables.push_back((unsigned int)vartype);
      }
      // In this case, also flag _where_ in the resulting list of variables the mapping starts
      m_beginMapped = dim;
    } else {
      int vartype = typeFromString(var);
      if (vartype < 0) {
	// Only flag the issue but otherwise take no action (assume non-argument use of a semicolon)
	std::cerr << "in CalibrationDataMappedHistogramContainer::computeVariableTypes(): cannot construct variable type from name "
		  << var << std::endl;
      } else {
	m_variables.push_back((unsigned int)vartype);
      }
    }
  }

  // After doing this, we should always have a non-null vector!
  assert(m_variables.size() > 0);

  // Also compute the validity bounds for this calibration object
  const_cast<CalibrationDataMappedHistogramContainer*>(this)->checkBounds();
}


// check the bounds of validity for this calibration object

void
CalibrationDataMappedHistogramContainer::checkBounds()
{
  const TH1* hist = dynamic_cast<const TH1*>(m_objResult);
  if (!hist) {
    std::cerr << "in CalibrationDataHistogramContainer::checkBounds(): object type does not derive from TH1" << std::endl;
    return;
  } else if (hist->GetDimension() + int(m_mapped.size()) - 1 != int(m_variables.size())) {
    std::cerr << "in CalibrationDataMappedHistogramContainer::checkBounds(): given number of variable types ("
	      << m_variables.size() << ") doesn't match (mapped) histogram dimension ("
	      << hist->GetDimension() + m_mapped.size() - 1 << ")!" << std::endl;
    return;
  }

  // Carry out the only cross-check that's possible for the binning: check that the dimensionality
  // for all bins matches the number of variables specified for the mapping
  for (unsigned int bin = 0; bin < m_bins.size(); ++bin)
    assert(m_bins[bin].getDimension() == m_mapped.size());

  for (unsigned int t = 0, t2 = 0; int(t) < hist->GetDimension(); ++t) {
    TAxis* axis;
    switch (t) {
    case 0: axis = hist->GetXaxis(); break;
    case 1: axis = hist->GetYaxis(); break;
    default: axis = hist->GetZaxis();
    }

    // Special case for the mapped dimension: here the only thing that can be done is to
    // cycle through all Bins and inspect the boundaries of each bin manually.
    if (t == m_beginMapped) {
      for (unsigned int mapped = 0; mapped < m_mapped.size(); ++mapped) {
	for (unsigned int bin = 0; bin < m_bins.size(); ++bin) {
	  double amax = m_bins[bin].getUpperBound(mapped), amin = m_bins[bin].getLowerBound(mapped);
	  if (bin == 0 || amax > m_upperBounds[m_variables[t2]]) m_upperBounds[m_variables[t2]] = amax;
	  if (bin == 0 || amin < m_lowerBounds[m_variables[t2]]) m_lowerBounds[m_variables[t2]] = amin;
	}
	++t2;
      }
    } else {
      // for (unsigned int t = 0; t < m_variables.size(); ++t) {
      //   if (m_variables[t] > m_upperBounds.size()) {
      // 	std::cerr << "in CalibrationDataHistogramContainer::checkBounds(): variable " << t << "type ("
      // 		  << m_variables[t] << "exceeds maximum type number (" << m_upperBounds.size() << ")!"
      // 		  << std::endl;
      // 	return;
      //   }
      // }
      double amax = axis->GetXmax(), amin = axis->GetXmin();
      if (amax < m_upperBounds[m_variables[t2]]) m_upperBounds[m_variables[t2]] = amax;
      if (amin > m_lowerBounds[m_variables[t2]]) m_lowerBounds[m_variables[t2]] = amin;
      ++t2;
    }
  }
}


// result retrieval

CalibrationDataContainer::CalibrationStatus
CalibrationDataMappedHistogramContainer::getResult(const CalibrationDataVariables& x,
						   double& result, TObject* obj) const
{
  if (!obj) {
    if (! m_objResult) {
      // std::cout << "retrieving central value pointer" << std::endl;
      m_objResult = GetValue("result");
    }
    obj = m_objResult;
  }
  TH1* hist = dynamic_cast<TH1*>(obj);
  if (! hist) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);
  // find the relevant "global" bin and retrieve its contents
  result = hist->GetBinContent(findBin());

  return inRange ? kSuccess : kRange;
}

// statistical uncertainty retrieval (special since it is stored with the result itself)

CalibrationDataContainer::CalibrationStatus
CalibrationDataMappedHistogramContainer::getStatUncertainty(const CalibrationDataVariables& x,
							    double& result) const
{
  if (! m_objResult) {
    // std::cout << "retrieving central value pointer" << std::endl;
    m_objResult = GetValue("result");
  }
  TH1* hist = dynamic_cast<TH1*>(m_objResult);
  if (! hist) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);
  // find the relevant "global" bin and retrieve its contents
  result = hist->GetBinError(findBin());

  return inRange ? kSuccess : kRange;
}

// general uncertainty retrieval

CalibrationDataContainer::CalibrationStatus
CalibrationDataMappedHistogramContainer::getUncertainty(const std::string& unc,
							const CalibrationDataVariables& x,
							UncertaintyResult& result, TObject* obj) const
{
  // treat statistical uncertainties separately (they are stored with the actual result)
  if (unc == "statistics") {
    double res;
    CalibrationStatus code = getStatUncertainty(x, res);
    if (code == kError) return code;
    result.first =   res;
    result.second = -res;
    return code;
  }

  if (!obj) obj = GetValue(unc.c_str());
  TH1* hist = dynamic_cast<TH1*>(obj);
  if (! hist) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);
  // find the relevant "global" bin and retrieve its contents
  Int_t bin = findBin();
  result.first = hist->GetBinError(bin);
  result.second = hist->GetBinError(bin);

  return inRange ? kSuccess : kRange;
}

/** test whether this calibration object is one for "continuous" calibration
    (this has some subtle consequences for the treatment of bin-to-bin correlations).
    The return value will be -1 in case this is not a "continuous" calibration object,
    and the axis number (0 for X, 1 for Y, 2 for Z) otherwise. */
int
CalibrationDataMappedHistogramContainer::getTagWeightAxis() const
{
  for (unsigned int type = 0; type < m_variables.size(); ++type)
    if (m_variables[type] == kTagWeight) {
      int hist_type = int(type);
      return (hist_type > int(m_beginMapped)) ? hist_type - m_mapped.size() + 1 : hist_type;
    }
  return -1;
}

// Set (by hand) the variables that will be mapped onto a single histogram axis

void
CalibrationDataMappedHistogramContainer::setMappedVariables(const std::vector<std::string>& variables)
{
  m_mapped = variables;
}

// List which variables get mapped onto a single histogram axis

std::vector<std::string>
CalibrationDataMappedHistogramContainer::getMappedVariables() const
{
  return m_mapped;
}

// Add bin to the present list
// Note the absence of a -1 in the return value: this is because ROOT's histogram axes start counting from 1

unsigned int
CalibrationDataMappedHistogramContainer::addBin(const Bin& bin)
{
  m_bins.push_back(bin);
  return m_bins.size();
}

// Return the number of mapped bins
// Note the absence of a -1 in the return value: this is because ROOT's histogram axes start counting from 1

unsigned int
CalibrationDataMappedHistogramContainer::getNMappedBins() const
{
  return m_bins.size();
}

// Find the mapped bin corresponding to the variables used for the mapping

Int_t
CalibrationDataMappedHistogramContainer::findMappedBin(double* x) const
{
  // First check quickly whether the last bin (cached) matches
  if (m_bins[m_lastBin].contains(x)) return m_lastBin + 1;

  // Search the whole array for a match
  for (unsigned int bin = 0; bin < m_bins.size(); ++bin)
    if (m_bins[bin].contains(x)) {
      m_lastBin = bin;
      return m_lastBin + 1;
    }
  std::cerr << "CalibrationDataMappedHistogramContainer::findMappedBin(): unable to find bin for mapping variables:";
  for (unsigned int d = 0; d < m_mapped.size(); ++d) std::cerr << "\t" << x[d];
  std::cerr << std::endl;
  // -1 means invalid..
  return -1;
}

// Find the bin corresponding to the computed variables (the computation is assumed to have just
// taken place and resulted in the m_vars array having been filled appropriately)

Int_t
CalibrationDataMappedHistogramContainer::findBin() const
{
  // Push the mapped variables onto an array.
  // Since we derive from TH1 this need never be more than 3 elements long.
  Int_t mapped[3];
  const TH1* hist = dynamic_cast<const TH1*>(m_objResult);
  Int_t ndim = hist->GetDimension();
  for (unsigned int dim = 0; dim < (unsigned int) ndim; ++dim) {
    if (dim == m_beginMapped) {
      if ((mapped[dim] = findMappedBin(&(m_vars[dim]))) < 0) return -1;
    } else {
      TAxis* axis;
      switch (dim) {
      case 0: axis = hist->GetXaxis(); break;
      case 1: axis = hist->GetYaxis(); break;
      default: axis = hist->GetZaxis();
      }
      mapped[dim] = axis->FindFixBin((dim < m_beginMapped) ? m_vars[dim] : m_vars[dim+m_mapped.size()-1]);
    }
  }
  return hist->GetBin(mapped[0], mapped[1], mapped[2]);
}

// Bin helper class methods

// default constructor (for persistency)
CalibrationDataMappedHistogramContainer::Bin::Bin():
  m_dimension(0), m_low(0), m_up(0)
{
}

// 'normal' constructor
CalibrationDataMappedHistogramContainer::Bin::Bin(unsigned int dimension, double* low, double* up):
  m_dimension(dimension)
{
  m_up  = new double[dimension];
  m_low = new double[dimension];
  for (unsigned int dim = 0; dim < dimension; ++dim) {
    m_up[dim]  = up[dim];
    m_low[dim] = low[dim];
  }
}

CalibrationDataMappedHistogramContainer::Bin::Bin(const CalibrationDataMappedHistogramContainer::Bin& other):
  m_dimension(other.m_dimension)
{
  m_up  = new double[m_dimension];
  m_low = new double[m_dimension];
  for (unsigned int dim = 0; dim < m_dimension; ++dim) {
    m_up[dim]  = other.m_up[dim];
    m_low[dim] = other.m_low[dim];
  }
}

CalibrationDataMappedHistogramContainer::Bin::~Bin()
{
  delete[] m_up;
  delete[] m_low;
}

bool
CalibrationDataMappedHistogramContainer::Bin::contains(double* x) const
{
  for (unsigned int dim = 0; dim < m_dimension; ++dim)
    if (x[dim] < m_low[dim] || x[dim] > m_up[dim]) return false;
  return true;
}

// ---------------------- CalibrationDataFunctionContainer class ----------------------

CalibrationDataFunctionContainer::CalibrationDataFunctionContainer(const char* name) :
  CalibrationDataContainer(name), m_objStatistics(0)
{
  // Reset the validity bounds to reflect 'no bounds'.

  m_lowerBounds.clear();
  m_lowerBounds.resize(maxParameters, -std::numeric_limits<double>::max());
  m_lowerBounds[kPt] = m_lowerBounds[kAbsEta] = 0;
  m_upperBounds.clear();
  m_upperBounds.resize(maxParameters, std::numeric_limits<double>::max());
}

CalibrationDataFunctionContainer::~CalibrationDataFunctionContainer()
{
}

// Determine which variable types are to be used for all objects (results + uncertainties).
// This needs to be done only once per calibration object, as the results will be cached
// (even if not persistified).

void
CalibrationDataFunctionContainer::computeVariableTypes() const
{
  if (! m_objResult) m_objResult = GetValue("result");

  std::string title(m_objResult->GetTitle());
  std::string::size_type pos = title.find(";");
  while (pos != std::string::npos && pos != title.size()) {
    title = title.substr(pos+1);
    pos = title.find(";");
    std::string var = title.substr(0, pos);
    int vartype = typeFromString(var);
    if (vartype < 0) {
      // Only flag the issue but otherwise take no action (assume non-argument use of a semicolon)
      std::cerr << "in CalibrationDataFunctionContainer::getVariableTypes(): cannot construct variable type from name "
		<< var << std::endl;
    } else {
      m_variables.push_back((unsigned int)vartype);
    }
  }

  // After doing this, we should always have a non-null vector!
  assert(m_variables.size() > 0);
}

// result retrieval

CalibrationDataContainer::CalibrationStatus
CalibrationDataFunctionContainer::getResult(const CalibrationDataVariables& x,
					    double& result, TObject* obj) const
{
  if (!obj) {
    if (! m_objResult) m_objResult = GetValue("result");
    obj = m_objResult;
  }
  TF1* func = dynamic_cast<TF1*>(obj);
  if (! func) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);
  result = func->EvalPar(m_vars);

  return inRange ? kSuccess : kRange;
}

// general uncertainty retrieval

CalibrationDataContainer::CalibrationStatus
CalibrationDataFunctionContainer::getUncertainty(const std::string& unc,
						 const CalibrationDataVariables& x,
						 UncertaintyResult& result, TObject* obj) const
{
  // treat statistical uncertainties separately (they are stored with the actual result)
  if (unc == "statistics") {
    double res;
    CalibrationStatus code = getStatUncertainty(x, res);
    if (code == kError) return code;
    result.first =   res;
    result.second = -res;
    return code;
  }

  if (!obj) obj = GetValue(unc.c_str());
  TF1* func = dynamic_cast<TF1*>(obj);
  if (! func) return kError;

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);

  // the "first" and "second" entries are filled with the
  // "positive" and "negative" uncertainties, respectively.
  // Note: no "negative" uncertainties implemented as yet!
  result.first = func->EvalPar(m_vars);
  result.second = -result.first;

  return inRange ? kSuccess : kRange;
}

// statistical uncertainty retrieval (special since the uncertainty is not parametrized)

CalibrationDataContainer::CalibrationStatus
CalibrationDataFunctionContainer::getStatUncertainty(const CalibrationDataVariables& x,
						     double& result) const
{
  // ensure that the requested objects exist
  if (! m_objResult) m_objResult = GetValue("result");
  TF1* func = dynamic_cast<TF1*>(m_objResult);
  if (! func) {
    // std::cerr << "... unable to retrieve the result" << std::endl;
    return kError;
  }

  if (! m_objStatistics) m_objStatistics = GetValue("statistics");
  //  m_objStatistics->Dump();
  TMatrixTSym<double>* cov = dynamic_cast<TMatrixTSym<double>*>(m_objStatistics);
  if (! cov) {
    // std::cerr << "... unable to retrieve the covariance matrix" << std::endl;
    return kError;
  }

  // select the relevant kinematic variables
  bool inRange = computeVariables(x);

  // use a large value for "eps": this multiplies the uncertainties that
  // are expected to be associated with the parameters. Choosing a large
  // value expresses the fact that we are not primarily interested in the
  // parabolic behaviour at the minimum
  // const Double_t eps = 1.0;
  // test: set to 0.5
  const Double_t eps = 0.5;

  int npar = func->GetNpar();
  if (npar == 0) {
    result = 0.;
    return kSuccess;
  }

  TMatrixT<double> gradients(npar,1);
  //  std::cout << "parametricVariance: gradients:";
  for (int ipar = 0; ipar < npar; ++ipar) {
    gradients(ipar,0) = func->GradientPar(ipar, m_vars, eps);
    // std::cout << " " << gradients(ipar,0);
  }
  //  std::cout << std::endl;

  // carry out the matrix multiplication
  TMatrixT<double> gradientsTransposed(TMatrixT<double>::kTransposed, gradients);
  //  std::cout << "parametricVariance: transposed gradients:";
  //  for (int ipar = 0; ipar < npar; ++ipar)
  //    std::cout << " " << gradients(0,ipar);
  //  std::cout << std::endl;
  TMatrixT<double> tmp1(*cov, TMatrixT<double>::kMult, gradients);
  //  std::cout << "parametricVariance: cov * gradients:";
  //  for (int ipar = 0; ipar < npar; ++ipar)
  //    std::cout << " " << tmp1(ipar,0);
  TMatrixT<double> tmp2(gradientsTransposed, TMatrixT<double>::kMult, tmp1);

  result = TMath::Sqrt(tmp2(0,0));

  return inRange ? kSuccess : kRange;
}
