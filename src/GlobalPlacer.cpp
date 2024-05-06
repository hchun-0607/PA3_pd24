#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

GlobalPlacer::GlobalPlacer(Placement &placement)
    : _placement(placement) {
}

void GlobalPlacer::place() {
    ////////////////////////////////////////////////////////////////////
    // This section is an example for analytical methods.
    // The objective is to minimize the following function:
    //      f(x,y) = 3*x^2 + 2*x*y + 2*y^2 + 7
    //
    // If you use other methods, you can skip and delete it directly.
    ////////////////////////////////////////////////////////////////////
    std::vector<Point2<double>> t(1);                   // Optimization variables (in this example, there is only one t)
    ExampleFunction foo(_placement);                    // Objective function
    const double kAlpha = 0.01;                         // Constant step size
    SimpleConjugateGradient optimizer(foo, t, kAlpha);  // Optimizer

    // Set initial point
    t[0] = 4.;  // This set both t[0].x and t[0].y to 4.

    // Initialize the optimizer
    optimizer.Initialize();

    // Perform optimization, the termination condition is that the number of iterations reaches 100
    // TODO: You may need to change the termination condition, which is determined by the overflow ratio.
    for (size_t i = 0; i < 50; ++i) {
        optimizer.Step();
        printf("iter = %3lu, f = %9.4f, x = %9.4f, y = %9.4f\n", i, foo(t), t[0].x, t[0].y);
    }

    ////////////////////////////////////////////////////////////////////
    // Global placement algorithm
    ////////////////////////////////////////////////////////////////////



    // TODO: Implement your global placement algorithm here.
    const size_t num_modules = _placement.numModules();  // You may modify this line.
    std::vector<Point2<double>> positions(num_modules);  // Optimization variables (positions of modules). You may modify this line.
    for(int i = 0; i < num_modules; i++){
        positions[i].x = rand()%1000;
        positions[i].y = rand()%1000;
    }
    bool check_wire = 0;
    bool check_density = 1;
    bool check_costfunction = 0;
    if(check_wire){
        Wirelength wirelength_(_placement);
        SimpleConjugateGradient optimizer_wire(wirelength_, positions, 0.01);
        optimizer_wire.Initialize();
        for(int i = 0; i < 50; i++){
            optimizer_wire.Step();
            printf("iter = %3lu, f = %9.4f\n, x = %9.4f, y = %9.4f", i, wirelength_(positions) , positions[0].x, positions[0].y);
        }
    }
    if(check_density){
        Density density_(_placement, 0.8, 100, 100);
        SimpleConjugateGradient optimizer_density(density_, positions, 0.01);
        optimizer_density.Initialize();
        //cout<<"check point1"<<endl;
        for(int i = 0; i < 50; i++){
            cout<<"i = " << i <<endl;
            optimizer_density.Step();
            printf("iter = %3lu, f = %9.4f\n, x = %9.4f, y = %9.4f\n", i,density_(positions) , positions[0].x, positions[0].y);
        }
    }
    if(check_costfunction){
        ObjectiveFunction cost_function(_placement, 0.1);
        SimpleConjugateGradient optimizer_cost(cost_function, positions, 0.01);
        optimizer_cost.Initialize();
        for(int i = 0; i < 50; i++){
            optimizer_cost.Step();
            cout<<i<<endl;
            printf("iter = %3lu, f = %9.4f\n, x = %9.4f, y = %9.4f", i,cost_function(positions) , positions[0].x, positions[0].y);
        } 
    }
    

    ////////////////////////////////////////////////////////////////////
    // Write the placement result into the database. (You may modify this part.)
    ////////////////////////////////////////////////////////////////////
    for (size_t i = 0; i < num_modules; i++) {
        _placement.module(i).setPosition(positions[i].x, positions[i].y);
    }
    cout<<"Global placement done!"<<endl;
}

void GlobalPlacer::plotPlacementResult(const string outfilename, bool isPrompt) {
    ofstream outfile(outfilename.c_str(), ios::out);
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl
            << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl
            << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT(outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop());
    outfile << "EOF" << endl;
    outfile << "# modules" << endl
            << "0.00, 0.00" << endl
            << endl;
    for (size_t i = 0; i < _placement.numModules(); ++i) {
        Module &module = _placement.module(i);
        plotBoxPLT(outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height());
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if (isPrompt) {
        char cmd[200];
        sprintf(cmd, "gnuplot %s", outfilename.c_str());
        if (!system(cmd)) {
            cout << "Fail to execute: \"" << cmd << "\"." << endl;
        }
    }
}

void GlobalPlacer::plotBoxPLT(ofstream &stream, double x1, double y1, double x2, double y2) {
    stream << x1 << ", " << y1 << endl
           << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl
           << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl
           << endl;
}
