
#include "IntegrationMethodTwoStep.h"
#include "NeighborList.h"

#ifndef __TWO_STEP_VICSEK_BASE__
#define __TWO_STEP_VICSEK_BASE__


class PYBIND11_EXPORT TwoStepVicsekBase : public IntegrationMethodTwoStep
    {
    public:
        //! Constructs the integration method and associates it with the system
        TwoStepVicsekBase(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<ParticleGroup> group,
                           std::shared_ptr<NeighborList> nlist,
                           Scalar r_cut,
                           Scalar v0,
                           Scalar delta,
                           Scalar eta,
                           unsigned int seed);
        virtual ~TwoStepVicsekBase();

    protected:
        std::shared_ptr<ParticleGroup> m_group;
        std::shared_ptr<NeighborList> m_nlist;
        Scalar m_r_cut;
        Scalar m_v0;
        Scalar m_delta;
        Scalar m_eta;
        unsigned int m_seed;
        Scalar m_deltaT;
    };

//! Exports the TwoStepLangevinBase class to python
void export_TwoStepVicsekBase(pybind11::module& m);

#endif // #ifndef __TWO_STEP_LANGEVIN_BASE__
