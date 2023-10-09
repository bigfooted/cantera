//! @file ReactorNet.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORNET_H
#define CT_REACTORNET_H

#include "Reactor.h"
#include "cantera/numerics/FuncEval.h"


namespace Cantera
{

class Array2D;
class Integrator;
class PreconditionerBase;

//! A class representing a network of connected reactors.
/*!
 *  This class is used to integrate the governing equations for a network of reactors
 *  that are time dependent (Reactor, ConstPressureReactor) connected by various
 *  means, for example Wall, MassFlowController, Valve, or PressureController; or
 *  reactors dependent on a single spatial variable (FlowReactor).
 *
 * @ingroup zerodGroup
 */
class ReactorNet : public FuncEval
{
public:
    ReactorNet();
    ~ReactorNet() override;
    ReactorNet(const ReactorNet&) = delete;
    ReactorNet& operator=(const ReactorNet&) = delete;

    //! @name Methods to set up a simulation
    //! @{

    //! Set the type of linear solver used in the integration.
    //! @param linSolverType type of linear solver. Default type: "DENSE"
    //! Other options include: "DIAG", "DENSE", "GMRES", "BAND"
    void setLinearSolverType(const string& linSolverType="DENSE");

    //! Set preconditioner used by the linear solver
    //! @param preconditioner preconditioner object used for the linear solver
    void setPreconditioner(shared_ptr<PreconditionerBase> preconditioner);

    //! Set the initial value of the independent variable (typically time).
    //! Default = 0.0 s. Restarts integration from this value using the current mixture
    //! state as the initial condition.
    void setInitialTime(CanteraDouble time);

    //! Get the initial value of the independent variable (typically time).
    /*!
     * @since New in %Cantera 3.0.
     */
    CanteraDouble getInitialTime() const {
        return m_initial_time;
    }

    //! Get the maximum integrator step.
    CanteraDouble maxTimeStep() const {
        return m_maxstep;
    }

    //! Set the maximum integrator step.
    void setMaxTimeStep(CanteraDouble maxstep);

    //! Set the maximum number of error test failures permitted by the CVODES
    //! integrator in a single step.
    void setMaxErrTestFails(int nmax);

    //! Set the relative and absolute tolerances for the integrator.
    void setTolerances(CanteraDouble rtol, CanteraDouble atol);

    //! Set the relative and absolute tolerances for integrating the
    //! sensitivity equations.
    void setSensitivityTolerances(CanteraDouble rtol, CanteraDouble atol);

    //! Current value of the simulation time [s], for reactor networks that are solved
    //! in the time domain.
    CanteraDouble time();

    //! Current position [m] along the length of the reactor network, for reactors that
    //! are solved as a function of space.
    CanteraDouble distance();

    //! Relative tolerance.
    CanteraDouble rtol() {
        return m_rtol;
    }

    //! Absolute integration tolerance
    CanteraDouble atol() {
        return m_atols;
    }

    //! Relative sensitivity tolerance
    CanteraDouble rtolSensitivity() const {
        return m_rtolsens;
    }

    //! Absolute sensitivity tolerance
    CanteraDouble atolSensitivity() const {
        return m_atolsens;
    }

    //! Problem type of integrator
    string linearSolverType() const;

    //! Returns the maximum number of internal integration steps the
    //! integrator will take before reaching the next output point
    int maxSteps();

    /**
     * Advance the state of all reactors in the independent variable (time or space).
     * Take as many internal steps as necessary to reach *t*.
     * @param t Time/distance to advance to (s or m).
     */
    void advance(CanteraDouble t);

    /**
     * Advance the state of all reactors in the independent variable (time or space).
     * Take as many internal steps as necessary towards *t*. If *applylimit* is true,
     * the advance step will be automatically reduced if needed to stay within limits
     * (set by setAdvanceLimit).
     * Returns the time/distance at the end of integration.
     * @param t Time/distance to advance to (s or m).
     * @param applylimit Limit advance step (boolean).
     */
    CanteraDouble advance(CanteraDouble t, bool applylimit);

    //! Advance the state of all reactors with respect to the independent variable
    //! (time or space). Returns the new value of the independent variable [s or m].
    CanteraDouble step();

    //! Add the reactor *r* to this reactor network.
    void addReactor(Reactor& r);

    //! Return a reference to the *n*-th reactor in this network. The reactor
    //! indices are determined by the order in which the reactors were added
    //! to the reactor network.
    Reactor& reactor(int n) {
        return *m_reactors[n];
    }

    //! Returns `true` if verbose logging output is enabled.
    bool verbose() const {
        return m_verbose;
    }

    //! Enable or disable verbose logging while setting up and integrating the
    //! reactor network.
    void setVerbose(bool v = true) {
        m_verbose = v;
        suppressErrors(!m_verbose);
    }

    //! Return a reference to the integrator. Only valid after adding at least one
    //! reactor to the network.
    Integrator& integrator();

    //! Update the state of all the reactors in the network to correspond to
    //! the values in the solution vector *y*.
    void updateState(CanteraDouble* y);

    //! Return the sensitivity of the *k*-th solution component with respect to
    //! the *p*-th sensitivity parameter.
    /*!
     *  The sensitivity coefficient @f$ S_{ki} @f$ of solution variable @f$ y_k
     *  @f$ with respect to sensitivity parameter @f$ p_i @f$ is defined as:
     *
     *  @f[ S_{ki} = \frac{1}{y_k} \frac{\partial y_k}{\partial p_i} @f]
     *
     *  For reaction sensitivities, the parameter is a multiplier on the forward
     *  rate constant (and implicitly on the reverse rate constant for
     *  reversible reactions) which has a nominal value of 1.0, and the
     *  sensitivity is nondimensional.
     *
     *  For species enthalpy sensitivities, the parameter is a perturbation to
     *  the molar enthalpy of formation, such that the dimensions of the
     *  sensitivity are kmol/J.
     */
    CanteraDouble sensitivity(size_t k, size_t p);

    //! Return the sensitivity of the component named *component* with respect to
    //! the *p*-th sensitivity parameter.
    //! @copydetails ReactorNet::sensitivity(size_t, size_t)
    CanteraDouble sensitivity(const string& component, size_t p, int reactor=0) {
        size_t k = globalComponentIndex(component, reactor);
        return sensitivity(k, p);
    }

    //! Evaluate the Jacobian matrix for the reactor network.
    /*!
     *  @param[in] t Time/distance at which to evaluate the Jacobian
     *  @param[in] y Global state vector at *t*
     *  @param[out] ydot Derivative of the state vector evaluated at *t*, with respect
     *      to *t*.
     *  @param[in] p sensitivity parameter vector (unused?)
     *  @param[out] j Jacobian matrix, size neq() by neq().
     */
    void evalJacobian(CanteraDouble t, CanteraDouble* y,
                      CanteraDouble* ydot, CanteraDouble* p, Array2D* j);

    // overloaded methods of class FuncEval
    size_t neq() const override {
        return m_nv;
    }

    size_t nReactors() const {
        return m_reactors.size();
    }

    void eval(CanteraDouble t, CanteraDouble* y, CanteraDouble* ydot, CanteraDouble* p) override;

    //! eval coupling for IDA / DAEs
    void evalDae(CanteraDouble t, CanteraDouble* y, CanteraDouble* ydot, CanteraDouble* p,
                 CanteraDouble* residual) override;

    void getState(CanteraDouble* y) override;
    void getStateDae(CanteraDouble* y, CanteraDouble* ydot) override;

    //! Return k-th derivative at the current state of the system
    virtual void getDerivative(int k, CanteraDouble* dky);

    void getConstraints(CanteraDouble* constraints) override;

    size_t nparams() const override {
        return m_sens_params.size();
    }

    //! Return the index corresponding to the component named *component* in the
    //! reactor with index *reactor* in the global state vector for the
    //! reactor network.
    size_t globalComponentIndex(const string& component, size_t reactor=0);

    //! Return the name of the i-th component of the global state vector. The
    //! name returned includes both the name of the reactor and the specific
    //! component, for example `'reactor1: CH4'`.
    string componentName(size_t i) const;

    //! Used by Reactor and Wall objects to register the addition of
    //! sensitivity parameters so that the ReactorNet can keep track of the
    //! order in which sensitivity parameters are added.
    //! @param name A name describing the parameter, for example the reaction string
    //! @param value The nominal value of the parameter
    //! @param scale A scaling factor to be applied to the sensitivity
    //!     coefficient
    //! @returns the index of this parameter in the vector of sensitivity
    //!     parameters (global across all reactors)
    size_t registerSensitivityParameter(const string& name, CanteraDouble value, CanteraDouble scale);

    //! The name of the p-th sensitivity parameter added to this ReactorNet.
    const string& sensitivityParameterName(size_t p) const {
        return m_paramNames.at(p);
    }

    //! Initialize the reactor network. Called automatically the first time
    //! advance or step is called.
    void initialize();

    //! Reinitialize the integrator. Used to solve a new problem (different
    //! initial conditions) but with the same configuration of the reactor
    //! network. Can be called manually, or automatically after calling
    //! setInitialTime or modifying a reactor's contents.
    void reinitialize();

    //! Called to trigger integrator reinitialization before further
    //! integration.
    void setNeedsReinit() {
        m_integrator_init = false;
    }

    //! Set the maximum number of internal integration steps the
    //! integrator will take before reaching the next output point
    //! @param nmax The maximum number of steps, setting this value
    //!             to zero disables this option.
    virtual void setMaxSteps(int nmax);

    //! Set absolute step size limits during advance
    void setAdvanceLimits(const CanteraDouble* limits);

    //! Check whether ReactorNet object uses advance limits
    bool hasAdvanceLimits() const;

    //! Retrieve absolute step size limits during advance
    bool getAdvanceLimits(CanteraDouble* limits) const;

    void preconditionerSetup(CanteraDouble t, CanteraDouble* y, CanteraDouble gamma) override;

    void preconditionerSolve(CanteraDouble* rhs, CanteraDouble* output) override;

    //! Get solver stats from integrator
    AnyMap solverStats() const;

    //! Set derivative settings of all reactors
    //! @param settings the settings map propagated to all reactors and kinetics objects
    virtual void setDerivativeSettings(AnyMap& settings);

protected:
    //! Check that preconditioning is supported by all reactors in the network
    virtual void checkPreconditionerSupported() const;

    void updatePreconditioner(CanteraDouble gamma) override;

    //! Estimate a future state based on current derivatives.
    //! The function is intended for internal use by ReactorNet::advance
    //! and deliberately not exposed in external interfaces.
    virtual void getEstimate(CanteraDouble time, int k, CanteraDouble* yest);

    //! Returns the order used for last solution step of the ODE integrator
    //! The function is intended for internal use by ReactorNet::advance
    //! and deliberately not exposed in external interfaces.
    virtual int lastOrder() const;

    vector<Reactor*> m_reactors;
    unique_ptr<Integrator> m_integ;

    //! The independent variable in the system. May be either time or space depending
    //! on the type of reactors in the network.
    CanteraDouble m_time = 0.0;

    //! The initial value of the independent variable in the system.
    CanteraDouble m_initial_time = 0.0;

    bool m_init = false;
    bool m_integrator_init = false; //!< True if integrator initialization is current
    size_t m_nv = 0;

    //! m_start[n] is the starting point in the state vector for reactor n
    vector<size_t> m_start;

    vector<CanteraDouble> m_atol;
    CanteraDouble m_rtol = 1.0e-9;
    CanteraDouble m_rtolsens = 1.0e-4;
    CanteraDouble m_atols = 1.0e-15;
    CanteraDouble m_atolsens = 1.0e-6;
    shared_ptr<PreconditionerBase> m_precon;
    string m_linearSolverType;

    //! Maximum integrator internal timestep. Default of 0.0 means infinity.
    CanteraDouble m_maxstep = 0.0;

    bool m_verbose = false;

    //! Indicates whether time or space is the independent variable
    bool m_timeIsIndependent = true;

    //! Names corresponding to each sensitivity parameter
    vector<string> m_paramNames;

    vector<CanteraDouble> m_ydot;
    vector<CanteraDouble> m_yest;
    vector<CanteraDouble> m_advancelimits;
    //! m_LHS is a vector representing the coefficients on the
    //! "left hand side" of each governing equation
    vector<CanteraDouble> m_LHS;
    vector<CanteraDouble> m_RHS;
};
}

#endif
