/**
 *  @file CVodesIntegrator.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CVODESWRAPPER_H
#define CT_CVODESWRAPPER_H

#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/SundialsContext.h"
#include "cantera/base/ctexceptions.h"

#include "sundials/sundials_nvector.h"

namespace Cantera
{

/**
 * Wrapper class for 'cvodes' integrator from LLNL.
 *
 * See FuncEval.h. Classes that use CVodesIntegrator:
 * ImplicitSurfChem, ReactorNet
 */
class CVodesIntegrator : public Integrator
{
public:
    /**
     *  Constructor. Default settings: dense Jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    CVodesIntegrator();
    ~CVodesIntegrator() override;
    void setTolerances(CanteraDouble reltol, size_t n, CanteraDouble* abstol) override;
    void setTolerances(CanteraDouble reltol, CanteraDouble abstol) override;
    void setSensitivityTolerances(CanteraDouble reltol, CanteraDouble abstol) override;
    void initialize(CanteraDouble t0, FuncEval& func) override;
    void reinitialize(CanteraDouble t0, FuncEval& func) override;
    void integrate(CanteraDouble tout) override;
    CanteraDouble step(CanteraDouble tout) override;
    CanteraDouble& solution(size_t k) override;
    CanteraDouble* solution() override;
    CanteraDouble* derivative(CanteraDouble tout, int n) override;
    int lastOrder() const override;
    int nEquations() const  override{
        return static_cast<int>(m_neq);
    }
    int nEvals() const override;
    void setMaxOrder(int n) override {
        m_maxord = n;
    }
    void setMethod(MethodType t) override;
    void setMaxStepSize(CanteraDouble hmax) override;
    void setMinStepSize(CanteraDouble hmin) override;
    void setMaxSteps(int nmax) override;
    int maxSteps() override;
    void setMaxErrTestFails(int n) override;
    AnyMap solverStats() const override;
    void setLinearSolverType(const string& linSolverType) override {
        m_type = linSolverType;
    }
    string linearSolverType() const override {
        return m_type;
    }
    void setBandwidth(int N_Upper, int N_Lower) override {
        m_mupper = N_Upper;
        m_mlower = N_Lower;
    }
    int nSensParams() override {
        return static_cast<int>(m_np);
    }
    CanteraDouble sensitivity(size_t k, size_t p) override;

    //! Returns a string listing the weighted error estimates associated
    //! with each solution component.
    //! This information can be used to identify which variables are
    //! responsible for integrator failures or unexpected small timesteps.
    string getErrorInfo(int N);

    //! Error message information provide by CVodes
    string m_error_message;

    CanteraDouble currentTime() const override {
        return m_time;
    }

    void setRootFunctionCount(size_t nroots) override;

protected:
    //! Applies user-specified options to the underlying CVODES solver. Called
    //! during integrator initialization or reinitialization.
    void applyOptions();

private:
    void sensInit(CanteraDouble t0, FuncEval& func);

    //! Check whether a CVODES method indicated an error. If so, throw an exception
    //! containing the method name and the error code stashed by the cvodes_err() function.
    void checkError(long flag, const string& ctMethod, const string& cvodesMethod) const;

    size_t m_neq = 0;
    void* m_cvode_mem = nullptr; //!< CVODES internal memory object
    SundialsContext m_sundials_ctx; //!< SUNContext object for Sundials>=6.0
    void* m_linsol = nullptr; //!< Sundials linear solver object
    void* m_linsol_matrix = nullptr; //!< matrix used by Sundials
    FuncEval* m_func = nullptr;
    CanteraDouble m_t0 = 0.0;

    //! The current system time, corresponding to #m_y
    CanteraDouble m_time;

    //! The latest time reached by the integrator. May be greater than #m_time.
    CanteraDouble m_tInteg;

    //! The system state at #m_time
    N_Vector m_y = nullptr;

    N_Vector m_abstol = nullptr;
    N_Vector m_dky = nullptr;
    string m_type = "DENSE";
    int m_itol;
    int m_method;
    int m_maxord = 0;
    CanteraDouble m_reltol = 1e-9;
    CanteraDouble m_abstols = 1e-15;
    CanteraDouble m_reltolsens = 1e-5;
    CanteraDouble m_abstolsens = 1e-4;
    size_t m_nabs = 0;
    CanteraDouble m_hmax = 0.0;
    CanteraDouble m_hmin = 0.0;
    int m_maxsteps = 20000;
    int m_maxErrTestFails = 0;
    N_Vector* m_yS = nullptr;
    size_t m_np = 0;
    int m_mupper = 0;
    int m_mlower = 0;
    //! Indicates whether the sensitivities stored in m_yS have been updated
    //! for at the current integrator time.
    bool m_sens_ok = false;
    size_t m_nRootFunctions = 0;
};

} // namespace

#endif
