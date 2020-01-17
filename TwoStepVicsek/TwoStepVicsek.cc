// Copyright (c) 2009-2019 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.


// Maintainer: joaander

//#include "hoomd/md/NeighborList.h"

#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
#include "hoomd/extern/pybind/include/pybind11/numpy.h"

#include <math.h>
#include <random>

#include "hoomd/RandomNumbers.h"
#include "hoomd/RNGIdentifiers.h"

#include "hoomd/VectorMath.h"
#include "hoomd/HOOMDMath.h"

#include "TwoStepVicsek.h"

#include <iostream>
namespace py = pybind11;
using namespace std;
using namespace hoomd;


TwoStepVicsek::TwoStepVicsek(
        std::shared_ptr<SystemDefinition> sysdef,
        std::shared_ptr<ParticleGroup> group,
        std::shared_ptr<NeighborList> nlist,
        Scalar r_cut,
        Scalar v0,
        Scalar delta,
        Scalar eta,
        unsigned int seed)
        : TwoStepVicsekBase(sysdef, group, nlist, r_cut, v0, delta, eta, seed),
          m_group(group), 
          m_nlist(nlist), 
          m_r_cut(r_cut),
          m_v0(v0), 
          m_delta(delta),
          m_eta(eta),
          m_seed(seed)
        {
        m_exec_conf->msg->notice(5) << "Constructing IntegrateVicsek" << endl;
        }

TwoStepVicsek::~TwoStepVicsek()
    {
    m_exec_conf->msg->notice(5) << "Destroying IntegrateVicsek" << endl;
    }


std::vector< std::string > TwoStepVicsek::getProvidedLogQuantities()
    {
    vector<string> result;
    return result;
    }

/*! \param quantity Name of the log quantity to get
    \param timestep Current time step of the simulation
    \param my_quantity_flag passed as false, changed to true if quantity logged here
*/
Scalar TwoStepVicsek::getLogValue(const std::string& quantity, unsigned int timestep, bool &my_quantity_flag)
    {
    return Scalar(0);
    }



void TwoStepVicsek::integrateStepOne(unsigned int timestep)
    {
    assert(m_nlist);
    // start by updating the neighborlist (see md/NeighborList.h for documentation)
    m_nlist->compute(timestep);
    ArrayHandle<unsigned int> h_n_neigh(m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(m_nlist->getHeadList(), access_location::host, access_mode::read);
    // retrieve position and orientation data for the system
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_orientation(m_pdata->getOrientationArray(), access_location::host, access_mode::readwrite);
    // calculate the aligned orientations of each particle
    const BoxDim& box = m_pdata->getGlobalBox();
    unsigned int group_size = m_group->getNumMembers();
    
    float orientations[group_size];
    for (unsigned int i = 0; i < group_size; i++)
        {
        // retreiving theta=Re[quaternion] (quat = x+yi+zj+wk) and converting to xy-cartesian
        float theta_i = h_orientation.data[i].x;
        float R_x = slow::cos(theta_i);
        float R_y = slow::sin(theta_i);
        Scalar3 r_i = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, 0);
        unsigned int N = 1;
        // loop over neighbors (O(N) neighbors list method)
        const unsigned int myHead = h_head_list.data[i];
        const unsigned int size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int k = 0; k < size; k++)
            {
            unsigned int j = h_nlist.data[myHead + k];
            Scalar3 r_j = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, 0);
            Scalar3 dr = r_i - r_j;
            // apply periodic boundary conditions
            dr = box.minImage(dr);
            float drsq = dot(dr,dr);
            float r_cut_sq = m_r_cut*m_r_cut;
            // if within cut-off radius, particle 'j' is a neighbor
            if(drsq < r_cut_sq && j != i)
                {
                float theta_j = h_orientation.data[j].x;
                R_x += slow::cos(theta_j);
                R_y += slow::sin(theta_j);
                N++;
                }
            }
        // Initialize the RNG
        ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
        unsigned int ptag = h_tag.data[i];
        RandomGenerator rng(RNGIdentifier::TwoStepLangevin, m_seed, ptag, timestep);
        // Create a uniform distribution on the unit circle
        hoomd::UniformDistribution<Scalar> uniform(Scalar(-3.14), Scalar(3.14));
        // add vectorial noise: N_i*eta*e^{i*phi}
        if(m_eta > 0)
            {
            Scalar phi = uniform(rng);
            R_x += N*m_eta*slow::cos(phi);
            R_y += N*m_eta*slow::sin(phi);
            }
        // align with neighbors
        float aligned_theta = atan2(R_y,R_x);
        // add scalar noise
        if(m_delta > 0)
            {
            Scalar varphi = uniform(rng);
            aligned_theta += m_delta*varphi;
            }
        // save the new orientation of particle 'i'
        orientations[i] = aligned_theta;
        }
    // update orientation of every particle (necessary as separate loop so that orientations are updated instantaneously)
    for (unsigned int i=0; i<group_size; i++)
        {
        h_orientation.data[i].x = orientations[i];
        }
    }


void TwoStepVicsek::integrateStepTwo(unsigned int timestep)
    {
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_orientation(m_pdata->getOrientationArray(), access_location::host, access_mode::readwrite);
    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::readwrite);
    const BoxDim& box = m_pdata->getBox();
    unsigned int group_size = m_group->getNumMembers();
    // for each particle in the system
    for (unsigned int i = 0; i < group_size; i++)
        {
        float thetai = h_orientation.data[i].x;
        // update the position of partice i
        h_pos.data[i].x += m_v0*slow::cos(thetai);
        h_pos.data[i].y += m_v0*slow::sin(thetai);
        h_pos.data[i].z = 0;
        // enforce period boundary condition
        box.wrap(h_pos.data[i], h_image.data[i]);
        }
    }


void export_TwoStepVicsek(py::module& m)
    {
    py::class_<TwoStepVicsek, std::shared_ptr<TwoStepVicsek> >(m, "TwoStepVicsek", py::base<TwoStepVicsekBase>())
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
