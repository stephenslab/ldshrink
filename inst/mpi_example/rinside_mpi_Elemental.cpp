#include <../include/LDshrink.h>
#include <El.hpp>
#include <mpi.h>     // mpi header
#include <RInside.h>// for the embedded R via RInside

#include <utility>
#include <chrono>
#include <thread>
using namespace El;

//[[Rcpp::depends(RcppEigen)]]

typedef Eigen::Map<Eigen::MatrixXd> mapmat;


int
main(int argc, char *argv[])
{
// mpi related
  El::Environment env( argc, argv );
  El::mpi::Comm comm = El::mpi::COMM_WORLD;
  RInside R(argc,argv);                      // create an embedded R instance
  
  El::Grid grid(comm);
  
  std::stringstream txt;
  R["m"]=85;
  R["Ne"]=1490.672741;
  R["cutoff"]=1e-3;
  R.parseEval("library(LDshrink); data(haplomat); data(mapdat)");
  mapmat hmat = Rcpp::as<mapmat>(R.parseEval("haplomat"));
  arrayxd_external map = Rcpp::as<arrayxd_external>(R.parseEval("mapdat"));
  arrayxd_external ldp = Rcpp::as<arrayxd_external>(R.parseEval("ldp <- c(m,Ne,cutoff,calc_theta(m))"));
  
  //  const int p = map.size();
  const int p=10;
  const int N=hmat.rows();
  //  const int csize =500;
  

  
  
  Eigen::ArrayXi chunka(3);
  
  const int local_row= grid.Row();
  const int local_col=grid.Col();

  
  const El::Int gridHeight = grid.Height();
  const El::Int gridWidth = grid.Width();
  const El::Int local_chunksize_row = El::Length( p, local_row, gridHeight );
  const El::Int local_chunksize_col = El::Length( p, local_col, gridWidth );
  // const El::Int local_chunksize_row = LDmat.LocalHeight();
  // const El::Int local_chunksize_col = LDmat.LocalWidth();

  
  const std::pair<int,int> Ablock= std::make_pair(local_row,local_chunksize_row);
  const std::pair<int,int> Bblock= std::make_pair(local_col,local_chunksize_col);
  

  const bool isDiag=local_row==local_col;


  const c_Matrix_external hmata(hmat.block(0,Ablock.first,N,Ablock.second).data(),N,Ablock.second);
  const c_Matrix_external hmatb(hmat.block(0,Bblock.first,N,Bblock.second).data(),N,Bblock.second);
  
  const c_arrayxd_external mapa(map.segment(Ablock.first,Ablock.second).data(),Ablock.second);
  const c_arrayxd_external mapb(map.segment(Bblock.first,Bblock.second).data(),Bblock.second);
//  El::Print(LDmat,"About to start calcLD_d\n");

  std::cout<<"I am of row: "<<local_row<<" col:"<<local_col<<" rowsize: "<<local_chunksize_row<<" colsize: "<<local_chunksize_col<<std::endl;
  Eigen::MatrixXd retmat =calcLD_d(hmata,mapa,hmatb,mapb,ldp,isDiag);    
  std::cout<<retmat<<std::endl;

  El::DistMatrix<double> LDmat(grid);
  LDmat.Attach(p,p,grid,0,0,retmat.data(),local_chunksize_col);
  El::Print(LDmat,"LDmat");
  El::Write(LDmat,"LDmat",El::ASCII);
//  El::Matrix<double> tretmat(retmat.rows(),retmat.cols(),&retmat.data());
  // 
  // 
  // 
  // R["txt"] = txt.str();	                // assign string var to R variable 'txt'
  // 
  // R.parseEvalQ("cat(txt)");                   // eval init string, ignoring any returns
  // 
  // MPI_Finalize();                             // mpi finalization
  // 
  // exit(0);
}
