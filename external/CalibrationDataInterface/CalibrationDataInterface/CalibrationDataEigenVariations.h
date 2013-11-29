//////////////////////////////////////////////////////////////////////
// CalibrationDataEigenVariations.h, (c) ATLAS Detector software
//////////////////////////////////////////////////////////////////////

#ifndef ANALYSISCALIBRATIONDATAINTERFACEEVVARIATIONS_H
#define ANALYSISCALIBRATIONDATAINTERFACEEVVARIATIONS_H

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "TMatrixDSym.h"

class TH1;

namespace Analysis 
{
  class CalibrationDataHistogramContainer;

  class CalibrationDataEigenVariations {
  public:
    CalibrationDataEigenVariations(const CalibrationDataHistogramContainer* cnt);
    ~CalibrationDataEigenVariations();

    /** exclude the source of uncertainty indicated by  name  from eigenvector calculations */
    void excludeNamedUncertainty(const std::string& name);
    /** carry out the eigenvector computations. Exclusion of any source of uncertainty has to
	be done before calling this method */
    void initialize();

    /** retrieve the number of named variations */
    unsigned int getNumberOfNamedVariations() const;
    /** list the named variations */
    std::vector<std::string> listNamedVariations() const;
    /** retrieve the integer index corresponding to the named variation.
	This can be used to speed up computations by avoiding string comparisons.
        Note that this function is not protected against passing an invalid name. */
    unsigned int getNamedVariationIndex(const std::string& name) const;

    /** retrieve the number of eigenvector variations */
    unsigned int getNumberOfEigenVariations() const;

    /** obtain the "up" and "down" variations for the given eigenvector number.
	The return value will be false if the eigenvector number is invalid. */
    bool getEigenvectorVariation(unsigned int variation, TH1*& up, TH1*& down) const;

    /** obtain the "up" and "down" variations for the named uncertainty.
	The return value will be false if the given name is not listed as
	being excluded from the eigenvector calculations. */
    bool getNamedVariation(const std::string& name, TH1*& up, TH1*& down) const;
    /** obtain the "up" and "down" variations for the source uncertainty
	pointed to by the given index (which is assumed to correspond to the
	value retrieved using getNamedVariationIndex()).
	The return value will be false if the index is out of bounds. */
    bool getNamedVariation(unsigned int nameIndex, TH1*& up, TH1*& down) const;

    /** also provide (some) access to the underlying information:
	covariance matrix corresponding to eigenvector variations */
    TMatrixDSym getEigenCovarianceMatrix() const;
    
  private:
    /** container object containing the basic information */
    const CalibrationDataHistogramContainer* m_cnt;

    /** flag whether the initialization has been carried out */
    bool m_initialized;

    /** named variations */
    std::map<std::string, unsigned int> m_namedIndices;
    std::vector<std::pair<TH1*, TH1*> > m_named;

    /** eigenvector variations */
    std::vector<std::pair<TH1*, TH1*> > m_eigen;

    // /** @ data members needed for eigenvector method **/
    // /** the map stores the int which is needed to access the other vector<> objects **/
    // mutable std::map<std::string, unsigned int> m_eigenvectorMethod_index;
    // std::vector<TMatrixT<double> > m_eigenvectorMethod_matrix;
    // std::vector<std::vector<TObject*> > m_eigenvectorMethod_uncUpProvider;
    // std::vector<std::vector<TObject*> > m_eigenvectorMethod_uncDownProvider;
    // std::vector<std::vector<TObject*> > m_eigenvectorMethod_uncProvider;      
  };

}

#endif // ANALYSISCALIBRATIONDATAINTERFACEEVVARIATIONS_H
