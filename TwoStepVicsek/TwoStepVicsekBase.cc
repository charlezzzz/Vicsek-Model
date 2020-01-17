
#include "TwoStepVicsekBase.h"

namespace py = pybind11;
using namespace std;

TwoStepVicsekBase::TwoStepVicsekBase(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<ParticleGroup> group,
                           std::shared_ptr<NeighborList> nlist,
                           Scalar r_cut,
                           Scalar v0,
                           Scalar delta,
                           Scalar eta,
                           unsigned int seed)
    : IntegrationMethodTwoStep(sysdef, group),
      m_group(group), 
      m_nlist(nlist), 
      m_r_cut(r_cut),
      m_v0(v0), 
      m_delta(delta),
      m_eta(eta),
      m_seed(seed)
    {
    }

TwoStepVicsekBase::~TwoStepVicsekBase()
    {
    }

void export_TwoStepVicsekBase(py::module& m)
    {
    py::class_<TwoStepVicsekBase, std::shared_ptr<TwoStepVicsekBase> >(m, "TwoStepVicsekBase", py::base<IntegrationMethodTwoStep>())
        .def(py::init< std::shared_ptr<SystemDefinition>,
                       std::shared_ptr<ParticleGroup>,
                       std::shared_ptr<NeighborList>,
                       Scalar,
                       Scalar,
                       Scalar,
                       Scalar,
                       unsigned int>())
        ;
    }

//void export_TwoStepVicsekBase(pybind11::module& m);

