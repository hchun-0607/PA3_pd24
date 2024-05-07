#define _GLIBCXX_USE_CXX11_ABI 0  // Align the ABI version to avoid compatibility issues with `Placment.h`
#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <vector>

#include "Placement.h"
#include "Point.h"

/**
 * @brief Base class for objective functions
 */
class BaseFunction {
   public:
    /////////////////////////////////
    // Conssutructors
    /////////////////////////////////

    BaseFunction(const size_t &input_size) : grad_(input_size) {}

    /////////////////////////////////
    // Accessors
    /////////////////////////////////

    const std::vector<Point2<double>> &grad() const { return grad_; }
    const double &value() const { return value_; }

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    // Forward pass, compute the value of the function
    virtual const double &operator()(const std::vector<Point2<double>> &input) = 0;

    // Backward pass, compute the gradient of the function
    virtual const std::vector<Point2<double>> &Backward() = 0;

   protected:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> grad_;  // Gradient of the function
    double value_;                      // Value of the function
};

/**
 * @brief Example function for optimization
 *
 * This is a simple example function for optimization. The function is defined as:
 *      f(t) = 3*t.x^2 + 2*t.x*t.y + 2*t.y^2 + 7
 */
class ExampleFunction : public BaseFunction {
   public:
    /////////////////////////////////
    // Constructors
    /////////////////////////////////

    ExampleFunction(Placement &placement);

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

   private:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> input_;  // Cache the input for backward pass
    Placement &placement_;
};

/**
 * @brief Wirelength function
 */
class Wirelength : public BaseFunction {
    // TODO: Implement the wirelength function, add necessary data members for caching
    

   public:
    /////////////////////////////////
    // Methods
    /////////////////////////////////
    Wirelength(Placement &placement);
    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

    private:
    double gamma = 1000.;
    std::vector<Point2<double>> input_;
    vector<double> exp_xi_gamma;
    vector<double> exp_yi_gamma;
    Placement &placement;
    size_t input_size;
    vector<vector<int>> pin_netlist;

};

/**
 * @brief Density function
 */
class Density : public BaseFunction {
    // TODO: Implement the density function, add necessary data members for caching
    
   public:
    /////////////////////////////////
    // Methods
    /////////////////////////////////
    Density(Placement &placement, double max_density, int num_bins_x, int num_bins_y);
    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

    private:
    double max_density = 2;
    vector<vector<double>> density_map;
    Placement &placement;
    int num_bins_x;
    int num_bins_y;
    vector<Point2<double>> input_;
    double bin_size_x;
    double bin_size_y;
    vector<vector<pair<double,double>>> center;
    vector<vector<double>> density_diff_map;
};

/**
 * @brief Objective function for global placement
 */
class ObjectiveFunction : public BaseFunction {
    // TODO: Implement the objective function for global placement, add necessary data
    // members for caching
    //
    // Hint: The objetive function of global placement is as follows:
    //       f(t) = wirelength(t) + lambda * density(t),
    // where t is the positions of the modules, and lambda is the penalty weight.
    // You may need an interface to update the penalty weight (lambda) dynamically.
   public:
    /////////////////////////////////
    // Methods
    ObjectiveFunction(Placement &placement, double lambda);
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input) override;
    const std::vector<Point2<double>> &Backward() override;

   private:
    Wirelength wirelength_;
    Density density_;
    double lambda_;
    Placement &placement_;
    vector<Point2<double>> input_;
    vector<Point2<double>> grad_wirelength;
    vector<Point2<double>> grad_density;
    int iter_num ;
};

#endif  // OBJECTIVEFUNCTION_H
