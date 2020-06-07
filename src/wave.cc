
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include<dune/grid/onedgrid.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>


#include"driver.hh"

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

 // Pročitaj ulaznu datoteku
    Dune::ParameterTree input_data;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("wave.ini",input_data);
    ptreeparser.readOptions(argc,argv,input_data);

    int refinement    = input_data.get<int>("grid.refinement");           //broj profinjenja
    double dt            =  input_data.get<double>("fem.dt");             //vremenski korak
    double c          =  input_data.get<double>("problem.speedofsound");  // brzina zvuka
    double T_end      =  input_data.get<double>("problem.T");             // konačno vrijeme
    std::string name  = input_data.get<std::string>("output.filename");

    //pravokutna mreža u 2D
  /*  constexpr int dim = 2;  // dimenzija mreže
    using GridType = Dune::YaspGrid<dim>;
    Dune::FieldVector<GridType::ctype,dim> L(1.0);             // Duljina stranice
    std::array<int,dim> s={10,10};          // broj ćelija po stranici
    GridType grid(L, s);*/

    //1D mreža
    constexpr int dim = 1;  // dimenzija mreže
    using GridType = Dune::OneDGrid;
    typedef Dune::OneDGrid::ctype DF;
    int N = 10;
    DF a = 0, b = 3*M_PI;
    GridType grid (N, a, b);

    if(refinement > 0)
         grid.globalRefine(refinement);

    auto gv = grid.leafGridView();
    driver(gv, dt, c, T_end, name);

    return 0;
}
