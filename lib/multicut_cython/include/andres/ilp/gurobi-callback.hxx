#pragma once
#ifndef ANDRES_ILP_GUROBI_CALLBACK_HXX
#define ANDRES_ILP_GUROBI_CALLBACK_HXX

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "gurobi_c++.h"

namespace andres {
namespace ilp {

class GurobiCallback {
public:
    enum class PreSolver {PRE_SOLVER_AUTO, PRE_SOLVER_PRIMAL, PRE_SOLVER_DUAL, PRE_SOLVER_NONE};
    enum class LPSolver {LP_SOLVER_PRIMAL_SIMPLEX, LP_SOLVER_DUAL_SIMPLEX, LP_SOLVER_BARRIER, LP_SOLVER_SIFTING};
    enum class Focus {FOCUS_FEASIBILITY, FOCUS_OPTIMALITY, FOCUS_BESTBOUND, FOCUS_BALANCED};

    class Callback: public GRBCallback
    {
    public:
        Callback(GurobiCallback& gurobi) :
            gurobi_(gurobi)
        {}

        virtual void separateAndAddLazyConstraints() = 0;

        void callback()
        try
        {
            if (where == GRB_CB_MIP)
            {
                objectiveBest_ = getDoubleInfo(GRB_CB_MIP_OBJBST);
            
                if(objectiveBest_ == GRB_INFINITY)
                    objectiveBest_ = std::numeric_limits<double>::infinity();

                objectiveBound_ = getDoubleInfo(GRB_CB_MIP_OBJBND);
            }
            else if (where == GRB_CB_MIPSOL)
                separateAndAddLazyConstraints();
        }
        catch (GRBException const& e)
        {
            throw std::runtime_error(std::string("GurobiCallback error ") + std::to_string(e.getErrorCode()) + ": " + e.getMessage() + " while executing GurobiCallback callback");
        }

        double label(size_t variableIndex)
        {
            return getSolution(gurobi_.gurobiModel_->getVar(variableIndex));
        }

        double getRuntime()
        {
            return getDoubleInfo(GRB_CB_RUNTIME);
        }

        double objective()
        {
            return getDoubleInfo(GRB_CB_MIPSOL_OBJ);
        }

    protected:
        template<class VariableIndexIterator, class CoefficientIterator>
        void addLazyConstraint(
            VariableIndexIterator viBegin,
            VariableIndexIterator viEnd,
            CoefficientIterator coefficient,
            double lowerBound,
            double upperBound
        ) 
        {
            GRBLinExpr expression;

            for (; viBegin != viEnd; ++viBegin, ++coefficient)
                expression += (*coefficient) * gurobi_.gurobiVariables_[static_cast<size_t>(*viBegin)];

            if (lowerBound == upperBound)
                addLazy(expression == lowerBound);
            else
            {
                if (lowerBound != -std::numeric_limits<double>::infinity())
                    addLazy(lowerBound <= expression);

                if (upperBound != std::numeric_limits<double>::infinity())
                    addLazy(expression <= upperBound);
            }
        }

        double objectiveBest_ { std::numeric_limits<double>::infinity() };
        double objectiveBound_ { -std::numeric_limits<double>::infinity() };

    private:
        GurobiCallback& gurobi_;
    };

    GurobiCallback();
    ~GurobiCallback();
    void setTimeLimit(const size_t);
    void setNumberOfThreads(const size_t);
    void setAbsoluteGap(const double);
    void setRelativeGap(const double);
    void setFocus(const Focus);
    void setVerbosity(const bool);
    void setLPSolver(const LPSolver);
    void setPreSolver(const PreSolver, const int = -1);
    void addVariables(const size_t, const double*);
    template<class Iterator>
        void setStart(Iterator);
    template<class VariableIndexIterator, class CoefficientIterator>
        void addConstraint(VariableIndexIterator, VariableIndexIterator,
                           CoefficientIterator, const double, const double);
    void setCallback(Callback&);
    void optimize();

    double objective() const;
    double bound() const;
    double gap() const;
    double label(const size_t) const;
    size_t numberOfThreads() const;
    double absoluteGap() const;
    double relativeGap() const;

private:
    GRBEnv gurobiEnvironment_;
    GRBModel* gurobiModel_ { nullptr };
    GRBVar* gurobiVariables_ { nullptr };
    GRBLinExpr gurobiObjective_;
    size_t nVariables_ { 0 };
};

inline
GurobiCallback::GurobiCallback()
:   gurobiEnvironment_()
{
    gurobiModel_ = new GRBModel(gurobiEnvironment_);
    setVerbosity(false);

    // aggressively seach for disconnected sub-problems:
    gurobiModel_->getEnv().set(GRB_IntParam_Disconnected, 2);
}

inline
GurobiCallback::~GurobiCallback() {
    if (gurobiModel_ != nullptr)
        delete gurobiModel_;

    if (gurobiVariables_ != nullptr)
        delete[] gurobiVariables_;
}

inline
void GurobiCallback::setTimeLimit(
    const size_t numberOfSeconds
) {
    gurobiModel_->getEnv().set(GRB_DoubleParam_TimeLimit, numberOfSeconds);
}

inline
void GurobiCallback::setNumberOfThreads(
    const size_t numberOfThreads
) {
    gurobiModel_->getEnv().set(GRB_IntParam_Threads, numberOfThreads);
}

inline
void GurobiCallback::setAbsoluteGap(
    double gap
) {
    gurobiModel_->getEnv().set(GRB_DoubleParam_MIPGapAbs, gap);
}

inline
void GurobiCallback::setRelativeGap(
    double gap
) {
    gurobiModel_->getEnv().set(GRB_DoubleParam_MIPGap, gap);
}

inline
void GurobiCallback::setFocus(
    const Focus focus
) {
    switch(focus) {
    case Focus::FOCUS_FEASIBILITY:
        gurobiModel_->getEnv().set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
        break;
    case Focus::FOCUS_OPTIMALITY:
        gurobiModel_->getEnv().set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_OPTIMALITY);
        break;
    case Focus::FOCUS_BESTBOUND:
        gurobiModel_->getEnv().set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BESTBOUND);
        break;
    default:
        gurobiModel_->getEnv().set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BALANCED);
        break;
    }
}

inline
void GurobiCallback::setVerbosity(
    const bool verbosity
) {
    if(verbosity) {
        gurobiModel_->getEnv().set(GRB_IntParam_OutputFlag, 1);
    }
    else {
        gurobiModel_->getEnv().set(GRB_IntParam_OutputFlag, 0);
    }
}

inline
void GurobiCallback::setPreSolver(
    const PreSolver preSolver,
    const int passes
) {
    switch(preSolver) {
    case PreSolver::PRE_SOLVER_NONE:
        gurobiModel_->getEnv().set(GRB_IntParam_Presolve, 0);
        return;
    case PreSolver::PRE_SOLVER_AUTO:
        gurobiModel_->getEnv().set(GRB_IntParam_PreDual, -1);
        break;
    case PreSolver::PRE_SOLVER_PRIMAL:
        gurobiModel_->getEnv().set(GRB_IntParam_PreDual, 0);
        break;
    case PreSolver::PRE_SOLVER_DUAL:
        gurobiModel_->getEnv().set(GRB_IntParam_PreDual, 1);
        break;
    }
    gurobiModel_->getEnv().set(GRB_IntParam_PrePasses, passes);

    // crushing allows the solver to translate variable indices in cuts 
    // to variable indices of the pre-solved problem
    /*
    if(crush) {
        gurobiEnvironment_.set(GRB_IntParam_PreCrush, 1);
    }
    else {
        gurobiEnvironment_.set(GRB_IntParam_PreCrush, 0);
    }
    */
}

inline
void GurobiCallback::setLPSolver(
    const LPSolver lpSolver
) {
    switch(lpSolver) {
    case LPSolver::LP_SOLVER_PRIMAL_SIMPLEX:
        gurobiModel_->getEnv().set(GRB_IntParam_NodeMethod, 0);
        break;
    case LPSolver::LP_SOLVER_DUAL_SIMPLEX:
        gurobiModel_->getEnv().set(GRB_IntParam_NodeMethod, 1);
        break;
    case LPSolver::LP_SOLVER_BARRIER:
        gurobiModel_->getEnv().set(GRB_IntParam_NodeMethod, 2);
        break;
    case LPSolver::LP_SOLVER_SIFTING:
        gurobiModel_->getEnv().set(GRB_IntParam_NodeMethod, 1); // dual simplex
        gurobiModel_->getEnv().set(GRB_IntParam_SiftMethod, 1); // moderate, 2 = aggressive
        break;
    }
}

inline
void GurobiCallback::addVariables(
    const size_t numberOfVariables,
    const double* coefficients
) {
    nVariables_ += numberOfVariables;
    gurobiVariables_ = gurobiModel_->addVars(numberOfVariables, GRB_BINARY);
    gurobiModel_->update();
    gurobiObjective_.addTerms(coefficients, gurobiVariables_, numberOfVariables);
    gurobiModel_->setObjective(gurobiObjective_);
}

inline
void GurobiCallback::setCallback(
    Callback& callback
) {
    gurobiModel_->getEnv().set(GRB_IntParam_LazyConstraints, 1);
    gurobiModel_->setCallback(&callback);
}

inline
void GurobiCallback::optimize() {
    gurobiModel_->optimize();
}

inline
double GurobiCallback::objective() const {
    return gurobiModel_->get(GRB_DoubleAttr_ObjVal);
}

inline
double GurobiCallback::bound() const {
    return gurobiModel_->get(GRB_DoubleAttr_ObjBound);
}

inline
double GurobiCallback::gap() const {
    return (objective() - bound()) / (1.0 + std::fabs(objective()));
}

inline
double GurobiCallback::label(
    const size_t variableIndex
) const {
    return gurobiVariables_[variableIndex].get(GRB_DoubleAttr_X);
}

inline
size_t GurobiCallback::numberOfThreads() const {
    return gurobiModel_->getEnv().get(GRB_IntParam_Threads);
}

inline
double GurobiCallback::absoluteGap() const {
    return gurobiModel_->getEnv().get(GRB_DoubleParam_MIPGapAbs);
}

inline
double GurobiCallback::relativeGap() const {
    return gurobiModel_->getEnv().get(GRB_DoubleParam_MIPGap);
}

template<class VariableIndexIterator, class CoefficientIterator>
inline
void GurobiCallback::addConstraint(
    VariableIndexIterator viBegin,
    VariableIndexIterator viEnd,
    CoefficientIterator coefficient,
    double lowerBound,
    double upperBound
) {
    GRBLinExpr expression;
    for(; viBegin != viEnd; ++viBegin, ++coefficient) {
        expression += (*coefficient) * gurobiVariables_[static_cast<size_t>(*viBegin)];
    }
    if(lowerBound == upperBound) {
        GRBLinExpr exact(lowerBound);
        gurobiModel_->addConstr(expression, GRB_EQUAL, exact);
    }
    else {
        if(lowerBound != -std::numeric_limits<double>::infinity()) {
            GRBLinExpr lower(lowerBound);
            gurobiModel_->addConstr(expression, GRB_GREATER_EQUAL, lower);
        }
        if(upperBound != std::numeric_limits<double>::infinity()) {
            GRBLinExpr upper(upperBound);
            gurobiModel_->addConstr(expression,GRB_LESS_EQUAL, upper);
        }
    }
}

template<class Iterator>
inline
void GurobiCallback::setStart(
    Iterator valueIterator
)
{
    for(size_t j = 0; j < nVariables_; ++j, ++valueIterator) {
        gurobiVariables_[j].set(GRB_DoubleAttr_Start, static_cast<double>(*valueIterator));
    }
}

} // namespace ilp
} // namespace andres

#endif
