
#include <mpi.h>     // mpi header
#include <RInside.h>// for the embedded R via RInside
#include <../include/LDshrink.h>
#include <utility>


//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(mpi)]]
typedef Eigen::Map<Eigen::MatrixXd> mapmat;

int main(int argc, char *argv[]){
  // mpi related
  int world_rank, world_size;                       // node information
  MPI_Init(&argc,&argv);                      // mpi initialization
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);     // obtain current node rank
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);   // obtain total nodes running.
  
  RInside R(argc,argv);                      // create an embedded R instance
  
  
  
  std::stringstream txt;
  R["m"]=85;
  R["Ne"]=1490.672741;
  R["cutoff"]=1e-3;
  R.parseEval("library(LDshrink); data(haplomat); data(mapdat)");
  mapmat hmat = Rcpp::as<mapmat>(R.parseEval("haplomat"));
  arrayxd_external map = Rcpp::as<arrayxd_external>(R.parseEval("mapdat"));
  arrayxd_external ldp = Rcpp::as<arrayxd_external>(R.parseEval("ldp <- c(m,Ne,cutoff,calc_theta(m))"));
  
  const int p = map.size();
  const int N=hmat.rows();
  const int csize =500;
  
  const int trank=world_rank;
  
  
  
  
  //    txt << "Hello from node " << myrank         // node information
  //	<< " of " << nodesize << " nodes!" << std::endl;

  
  const std::vector<int> chunkrow={0,0,1,1};
  const std::vector<int> chunkcol={0,1,0,1};
  
  

    Eigen::ArrayXi chunka(3);
    if(world_rank>1){
      chunka(0)=chunkrow[world_rank-2];
      chunka(1)=chunkcol[world_rank-2];
      chunka(2)=csize;
    }else{
      chunka(0)=0;
      chunka(1)=0;
      chunka(2)=p;
    }
    
    arrayxi_external chunkmap(chunka.data(),chunka.size());
    const std::pair<int,int> Ablock=chunk_block(chunka(0),chunka(2),p);  
    const std::pair<int,int> Bblock=chunk_block(chunka(1),chunka(2),p);  
    txt << world_rank << " of "<< world_size <<" starts at " <<Ablock.first<<" and is of size "<<Ablock.second<<std::endl;    
    const c_Matrix_external hmata(hmat.block(0,Ablock.first,N,Ablock.second).data(),N,Ablock.second);
    const c_Matrix_external hmatb(hmat.block(0,Bblock.first,N,Bblock.second).data(),N,Bblock.second);
    const c_arrayxd_external mapa(map.segment(Ablock.first,Ablock.second).data(),Ablock.second);
    const c_arrayxd_external mapb(map.segment(Bblock.first,Bblock.second).data(),Bblock.second);
    Eigen::MatrixXd retmat =calcLD_d(hmata,mapa,hmatb,mapb,ldp,Ablock.first==Bblock.first);
  
  
  
  R["txt"] = txt.str();	                // assign string var to R variable 'txt'
  
  R.parseEvalQ("cat(txt)");                   // eval init string, ignoring any returns
  
  MPI_Finalize();                             // mpi finalization
  
  exit(0);
}
