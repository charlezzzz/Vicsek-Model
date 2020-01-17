

#include "NeighborList.h"
#include "TwoStepVicsekBase.h"

#ifndef __TWO_STEP_VICSEK_H__
#define __TWO_STEP_VICSEK_H__

#include <hoomd/extern/pybind/include/pybind11/pybind11.h>


class PYBIND11_EXPORT TwoStepVicsek : public TwoStepVicsekBase
    {
    public:
        //! Constructor
        TwoStepVicsek(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<ParticleGroup> group,
                           std::shared_ptr<NeighborList> nlist,
                           Scalar r_cut,
                           Scalar v0,
                           Scalar delta,
                           Scalar eta,
                           unsigned int seed);

        virtual ~TwoStepVicsek();

        //! Returns a list of log quantities this integrator calculates
        virtual std::vector< std::string > getProvidedLogQuantities();

        //! Returns logged values
        Scalar getLogValue(const std::string& quantity, unsigned int timestep, bool &my_quantity_flag);

        //! Performs the second step of the integration
        virtual void integrateStepOne(unsigned int timestep);

        //! Performs the second step of the integration
        virtual void integrateStepTwo(unsigned int timestep);


    protected:
      std::shared_ptr<ParticleGroup> m_group;
      std::shared_ptr<NeighborList> m_nlist;
      Scalar m_r_cut;
      Scalar m_v0;
      Scalar m_delta;
      Scalar m_eta;
      unsigned int m_seed;
      Scalar m_deltaT;  //!< The time step
    };

//! Exports the IntegrateVicsek class to python
void export_TwoStepVicsek(pybind11::module& m);

#endif // #ifndef __TWO_STEP_VICSEK_H__
