
#include "MyTMINLP2.hpp"
#include "BonAmplInterface.hpp"

MyTMINLP::MyTMINLP(int PathNum, const std::vector< std::vector<int> > & vvEMatrix, const std::vector< std::vector<double> > & vvSMatrix) : printSol_(false)
{
    iNumberPath      = PathNum;  // m in draft
    vvCurrentEMatrix = vvEMatrix; 
    vvCurrentSMatrix = vvSMatrix;
    iNumberVertice    = (int) vvCurrentEMatrix.size(); // n in draft
    if (iNumberVertice != (int) vvCurrentSMatrix.size())
    {
        cerr<<"sizes of matrix don't match!"<<endl;
        exit(1);
    }
    iNumberVariables = iNumberPath * iNumberVertice; // mn in draft
}

bool 
MyTMINLP::get_variables_types(Index n, VariableType* var_types)
{
    std::cout<<"get_variables_types n="<<n<<std::endl;
    
    int i;
//    for (i=0; i < iNumberVariables; i++)
//        var_types[i] = BINARY;
    var_types[0] = BINARY;
    return true;
  //var_types[0] = BINARY;
  //var_types[1] = CONTINUOUS;
  //var_types[2] = CONTINUOUS;
  //var_types[3] = INTEGER;
  //return true;
}


bool 
MyTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
    std::cout<<"get_variables_linearity n="<<n<<std::endl;
    int i;
//    for (i=0; i < iNumberVariables; i++)
//        var_types[i] = Ipopt::TNLP::NON_LINEAR;
    var_types[0] = Ipopt::TNLP::LINEAR; //change in the future
    return true;
  //var_types[0] = Ipopt::TNLP::LINEAR;
  //var_types[1] = Ipopt::TNLP::NON_LINEAR;
  //var_types[2] = Ipopt::TNLP::NON_LINEAR;
  //var_types[3] = Ipopt::TNLP::LINEAR;
  //return true;
}


bool 
MyTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
    std::cout<<"get_constraints_linearity m="<<m<<std::endl;
    const_types[0] = Ipopt::TNLP::LINEAR;

    const_types[1] = Ipopt::TNLP::LINEAR;
   // const_types[2] = Ipopt::TNLP::NON_LINEAR;

   // const_types[3] = Ipopt::TNLP::NON_LINEAR;
    const_types[2] = Ipopt::TNLP::LINEAR;

/*
    assert(m==iNumberVertice-2+2*iNumberVariables-2*iNumberPath);
    int i;
    for (i=0; i<(iNumberVertice-2); i++)
        const_types[i] = Ipopt::TNLP::LINEAR;
    for (i=iNumberVertice-2; i< iNumberVertice-2+2*iNumberVariables-2*iNumberPath; i++)
        const_types[i] = Ipopt::TNLP::NON_LINEAR;
*/
    return true;
  //assert (m==3);
  //const_types[0] = Ipopt::TNLP::NON_LINEAR;
  //const_types[1] = Ipopt::TNLP::LINEAR;
  //const_types[2] = Ipopt::TNLP::LINEAR;
  //return true;
}
bool 
MyTMINLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g,
                       Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
    std::cout<<"get_nlp_info"<<std::endl;
    n = 1;
    m = 3;
    nnz_jac_g = 3;
    nnz_h_lag = 0;
 /*   n = iNumberVariables;
    m = iNumberVertice-2+2*iNumberVariables-2*iNumberPath;
    nnz_jac_g = 2*iNumberVariables*iNumberVertice - iNumberVariables - 2*iNumberPath;
    nnz_h_lag = iNumberVariables*(iNumberVertice+1)/2;*/
    index_style = Ipopt::TNLP::C_STYLE;
    return true;
//  n = 4;//number of variable
//  m = 3;//number of constraints
//  nnz_jac_g = 7;//number of non zeroes in Jacobian
//  nnz_h_lag = 2;//number of non zeroes in Hessian of Lagrangean
//  index_style = TNLP::FORTRAN_STYLE;
//  return true;
}

bool 
MyTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
    std::cout<<"get_bounds_info n="<<n<<" m="<<m<<std::endl;
    x_l[0] = 0.;
    x_u[0] = 1.;

    g_l[0] = 1;
    g_u[0] = DBL_MAX;

    g_l[1] = g_u[1] = 0.;
    g_l[2] = g_u[2] = 0.;

//    g_l[3] = g_u[3] = 0.;
//    g_l[4] = g_u[4] = 0.;

/*    assert(n==iNumberVariables);
    assert(m==(2*iNumberVariables - 2*iNumberPath + iNumberVertice - 2));

    int i, j; 
    for (i=0; i<iNumberPath; i++)
    {
        x_l[iNumberVertice * i + 1 - 1] = 1;
        x_u[iNumberVertice * i + 1 - 1] = 1;
        x_l[iNumberVertice * i + iNumberVertice - 1] = 1;
        x_u[iNumberVertice * i + iNumberVertice - 1] = 1;

        for (j=2; j<iNumberVertice; j++)
        {
            x_l[iNumberVertice * i + j - 1] = 0;
            x_u[iNumberVertice * i + j - 1] = 1;
        }
    }

    for (i=0; i < (iNumberVertice-2); i++)
    {
        g_l[i] = 1;
        g_u[i] = DBL_MAX;
    }

    for (i=(iNumberVertice-2); i<(2*iNumberVariables - 2*iNumberPath + iNumberVertice - 2); i++)
    {
        g_l[i]=g_u[i]=0;
    }*/

    return true;

  //assert(n==4);
  //assert(m==3);
  //x_l[0] = 0.;
  //x_u[0] = 1.;
  
  //x_l[1] = 0.;
  //x_u[1] = DBL_MAX;
  
  //x_l[2] =0.;
  //x_u[2] = DBL_MAX;
  
  //x_l[3] = 0;
  //x_u[3] = 5;
  
  //g_l[0] = -DBL_MAX;
  //g_u[0] = 1./4.;

  //g_l[1] = -DBL_MAX;
  //g_u[1] = 0;
  
  //g_l[2] = -DBL_MAX;
  //g_u[2] = 2;
  //return true;
}

bool 
MyTMINLP::get_starting_point(Index n, bool init_x, Number* x,
                             bool init_z, Number* z_L, Number* z_U,
                             Index m, bool init_lambda,
                             Number* lambda)
{
    std::cout<<"get_starting_point n="<<n<<" init_x="<<init_x<<" init_z="<<init_z<<" m="<<m<<" init_lambda="<<init_lambda<<std::endl;
    x[0] = 1;

/*      assert(n==iNumberVariables);
    assert(m==(2*iNumberVariables - 2*iNumberPath + iNumberVertice - 2));
    assert(init_x);
    assert(!init_lambda);

    int i;
    for (i=0; i<iNumberVariables; i++)
        x[i] = 1;
*/
    return true;
//  assert(n==4);
//  assert(m==3);
  
//  assert(init_x);
//  assert(!init_lambda);
//  x[0] = 0;
//  x[1] = 0;
//  x[2] = 0;
//  x[3] = 0;
//  return true;
}

bool 
MyTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    std::cout<<"eval_f n="<<n<<" new_x="<<new_x<<"x_0="<<x[0]<<std::endl;

    int j, k;
    obj_value = 0;
    
    obj_value = obj_value - 1 * x[0] - 1 * x[0]  ;
  /*    assert(n==iNumberVariables);

    int i, j, k;
    obj_value = 0;
    for (i=0; i<iNumberPath; i++)
        for (j=0; j<iNumberVertice; j++)
            for (k=0; k<iNumberVertice; k++)
            {
                obj_value = obj_value - vvCurrentSMatrix[i][j] * x[i*iNumberVertice+j] * x[i*iNumberVertice+k];
            } */
    return true;
  //assert(n==4);
  //obj_value = - x[0] - x[1] - x[2];
  //return true;
}

bool
MyTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    std::cout<<"eval_grad_f n="<<n<<" new_x="<<new_x<<" x_0="<<x[0]<<std::endl;
    int i;
    for (i=0; i<1; i++)
        grad_f[i] = 0;
  //  for (i=0; i<3; i++)
  //  {
        grad_f[0] =0  - 1 - 1 ;
  //  }
/*
 assert(n==iNumberVariables);

    int i, j, k, iXIndex;

    for (i=0; i<iNumberPath; i++)
        for (j=0; j<iNumberVertice; j++)
        {
            iXIndex = i*iNumberVertice+j;
            grad_f[iXIndex] = 0;
            for (k=0; k<iNumberVertice; k++)
            {
                grad_f[iXIndex] = grad_f[iXIndex] - 2 * vvCurrentSMatrix[j][k] * x[i*iNumberVertice+k];
            }
        }
*/
    return true;
 // assert(n==4);
 // grad_f[0] = -1.;
 // grad_f[1] = -1.;  
 // grad_f[2] = -1.;
 // grad_f[3] = 0.;
 // return true;
}

bool
MyTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    std::cout<<"eval_g n="<<n<<" new_n="<<new_x<<" m="<<m<<" x_0="<<x[0]<<std::endl;

    int i;

    g[0] = x[0];
  
    g[1] = -1 +  x[0]  ;
    g[1] = -1 + x[0];
   // g[2] = -x[0] + vvCurrentEMatrix[1][2] * x[0]  + vvCurrentEMatrix[1][1] * x[0] * x[0] + vvCurrentEMatrix[1][0] * x[0] ;

  //  g[3] = -x[0] + vvCurrentEMatrix[0][1] * x[0] + vvCurrentEMatrix[1][1] * x[0] * x[0] + vvCurrentEMatrix[2][1] * x[0];
 //   g[4] = -1 + vvCurrentEMatrix[0][2] + vvCurrentEMatrix[1][2] * x[0] + vvCurrentEMatrix[2][2];

    /*
    assert(n==iNumberVariables);
    assert(m==(2*iNumberVariables - 2*iNumberPath + iNumberVertice - 2));  

    int i,j,k,iGIndex;
    for (i=0; i<(iNumberVertice-2); i++)
    {
        g[i] = 0;
        for(j=0; j<iNumberPath; j++)
            g[i] += x[j*iNumberVertice+i+1] ; // x[ij] in draft with 2<=j<=(n-1)
    }
    for (i=0; i<iNumberPath; i++)
        for (j=0; j<(iNumberVertice-1); j++)
        {
            iGIndex = iNumberVertice - 2 + i*(iNumberVertice-1) + j;
            g[iGIndex] = -x[i*iNumberVertice+j];
            for (k=0; k<iNumberVertice; k++)
                g[iGIndex] += x[i*iNumberVertice+j] * x[i*iNumberVertice+k] * vvCurrentEMatrix[j][k];
        }
    for (i=0; i<iNumberPath; i++)
        for (j=1; j<iNumberVertice; j++)
        {
            iGIndex = iNumberVertice - 2 + iNumberPath*(iNumberVertice-1) + i*(iNumberVertice-1)+j-1;
            g[iGIndex] = -x[i*iNumberVertice+j];
            for (k=0; k<iNumberVertice; k++)
                g[iGIndex] += x[i*iNumberVertice+j] * x[i*iNumberVertice+k] * vvCurrentEMatrix[k][j];
        }*/
    return true;
//  assert(n==4);
//  assert(m==3);
  
//  g[0] = (x[1] - 1./2.)*(x[1] - 1./2.) + (x[2] - 1./2.)*(x[2] - 1./2.);
//  g[1] = x[0] - x[1];
//  g[2] = x[0] + x[2] + x[3];
  
//  return true;
}

bool
MyTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
                     Index m, Index nnz_jac, Index* iRow, Index *jCol,
                     Number* values)
{
    if (values == NULL)
    {
        std::cout<<"eval_jac_g n="<<n<<" new_x="<<new_x<<" m="<<m<<" nnz_jac="<<nnz_jac<<std::endl;
        iRow[0] = 0;
        jCol[0] = 0;

        iRow[1] = 1;
        jCol[1] = 0;
        iRow[2] = 2;
        jCol[2] = 0;
   //     iRow[3] = 3;
   //     jCol[3] = 0;

   //     iRow[4] = 4;
   //     jCol[4] = 0;



    }else
    {
        std::cout<<"eval_jac_g n="<<n<<" new_x="<<new_x<<" m="<<m<<" nnz_jac="<<nnz_jac;
        std::cout<<" x_0="<<x[0]<<std::endl;
        values[0] = 1;
        values[1] = 1;
        values[2] = 1;
     //   values[2] = -1+vvCurrentEMatrix[1][2]+vvCurrentEMatrix[1][0]+2*vvCurrentEMatrix[1][1]*x[0];
     //   values[3] = -1 + vvCurrentEMatrix[0][1] + 2*vvCurrentEMatrix[1][1]*x[0] +  vvCurrentEMatrix[2][1] ;
     //   values[4] = vvCurrentEMatrix[1][2];

    }


    /*
    assert(n==iNumberVariables);
    assert(nnz_jac == (2*iNumberVariables*iNumberVertice - iNumberVariables - 2*iNumberPath));

    int i, j, k, iEleNumber=0;
    if (values == NULL)
    {
        iEleNumber = 0;
        for (j=1; j<=(iNumberVertice-2); j++)
            for (i=0; i<iNumberPath; i++)
            {
                iRow[iEleNumber] = j-1;
                jCol[iEleNumber] = i*iNumberVertice + j;
                iEleNumber++;
            }
        for (i=0; i<iNumberPath; i++)
            for (j=0; j<(iNumberVertice-1); j++)
                for (k=0; k<iNumberVertice; k++)
                {
                    iRow[iEleNumber] = iNumberVertice-2 + i*(iNumberVertice-1)+j;
                    jCol[iEleNumber] = i*iNumberVertice + k;
                    iEleNumber++;
                }
        for (i=0; i<iNumberPath; i++)
            for (j=1; j<iNumberVertice; j++)
                for (k=0; k<iNumberVertice; k++)
                {
                    iRow[iEleNumber] = iNumberVertice-2 + iNumberPath*(iNumberVertice-1) + i*(iNumberVertice-1)+j-1;
                    jCol[iEleNumber] = i*iNumberVertice + k;
                    iEleNumber++;
                }
        cout<<"eval_jac_g_1: "<<iEleNumber<<endl;
    } else
    {
        iEleNumber = 0;
        for (j=1; j<=(iNumberVertice-2); j++)
            for (i=0; i<iNumberPath; i++)
            {
                values[iEleNumber] = 1;
                iEleNumber++;
            }
        for (i=0; i<iNumberPath; i++)
            for (j=0; j<(iNumberVertice-1); j++)
                for (k=0; k<iNumberVertice; k++)
                {
                    if (k!=j)
                    {
                        values[iEleNumber] = vvCurrentEMatrix[j][k] * x[i*iNumberVertice + j];
                    }else
                    {
                        values[iEleNumber] = 2*vvCurrentEMatrix[j][j]*x[i*iNumberVertice + j] - 1;
                    }
                    iEleNumber++;
                }
        for (i=0; i<iNumberPath; i++)
            for (j=1; j<iNumberVertice; j++)
                for (k=0; k<iNumberVertice; k++)
                {
                    if (k!=j)
                    {
                        values[iEleNumber] = vvCurrentEMatrix[k][j] * x[i*iNumberVertice + j];
                    }else
                    {
                        values[iEleNumber] = 2*vvCurrentEMatrix[j][j] * x[i*iNumberVertice + j] - 1;
                    }
                    iEleNumber++;
                }
        cout<<"eval_jac_g_2: "<<iEleNumber<<endl;
    }*/
    return true;
/*  assert(n==4);
  assert(nnz_jac == 7);
  if(values == NULL) {
    iRow[0] = 2;
    jCol[0] = 1;

    iRow[1] = 3;
    jCol[1] = 1;
    
    iRow[2] = 1;
    jCol[2] = 2;
    
    iRow[3] = 2;
    jCol[3] = 2;
    
    iRow[4] = 1;
    jCol[4] = 3;
        
    iRow[5] = 3;
    jCol[5] = 3;
    
    iRow[6] = 3;
    jCol[6] = 4;
    return true;
  }
  else {
    values[0] = 1.;
    values[1] = 1;

    values[2] = 2*x[1] - 1;
    values[3] = -1.;

    values[4] = 2*x[2] - 1;
    values[5] = 1.;
    
    values[6] = 1.;
    
    return true;
  }*/
}

bool
MyTMINLP::eval_h(Index n, const Number* x, bool new_x,
                 Number obj_factor, Index m, const Number* lambda,
                 bool new_lambda, Index nele_hess, Index* iRow,
                 Index* jCol, Number* values)
{
    if (values == NULL)
    {
        std::cout<<"eval_h n="<<n<<" new_x="<<new_x<<" obj_factor="<<obj_factor;
        std::cout<<" m="<<m<<" new_lambda="<<new_lambda<<" nele_hess="<<nele_hess<<std::endl;
        //std::cout<<" x_1="<<x[1]<<std::endl;
  //      iRow[0] = 0;
  //      jCol[0] = 0;

    }else
    {
        std::cout<<"eval_h n="<<n<<" new_x="<<new_x<<" obj_factor="<<obj_factor;
        std::cout<<" m="<<m<<" new_lambda="<<new_lambda<<" nele_hess="<<nele_hess;
        std::cout<<" x_0="<<x[0];
        std::cout<<" l_0="<<lambda[0]<<" l_1="<<lambda[1]<<" l_2="<<lambda[2];
        std::cout<<" l_3="<<lambda[3]<<" l_4="<<lambda[4]<<std::endl;


    //    values[0] = obj_factor * (-2) * vvCurrentSMatrix[1][1];

    //    values[0] += lambda[2] * 2 * vvCurrentEMatrix[1][1];

    //    values[0] += lambda[3] * 2 * vvCurrentEMatrix[1][1];

    }

    /*
    assert(n==iNumberVariables);
    assert(m==(2*iNumberVariables - 2*iNumberPath + iNumberVertice - 2));
    assert(nele_hess==(iNumberVariables*(iNumberVertice+1)/2));

    int i, j, k, iEleNumber=0;

    if (values == NULL)
    {
        iEleNumber = 0;
        for (i=0; i<iNumberPath; i++)
            for (j=0; j<iNumberVertice; j++)
                for (k=0; k<=j; k++)
                {
                    iRow[iEleNumber] = i*iNumberVertice + j;
                    jCol[iEleNumber] = i*iNumberVertice + k;
                    iEleNumber++;
                }
        cout<<"eval_h_1: "<<iEleNumber<<endl;
    }else
    {
        iEleNumber = 0;
        for (i=0; i<iNumberPath; i++)
            for (j=0; j<iNumberVertice; j++)
                for (k=0; k<=j; k++)
                {
                    if (k!=j)
                    {
                        values[iEleNumber] = obj_factor * (-2) * vvCurrentSMatrix[j][k];
                        if (j<(iNumberVertice-1))
                            values[iEleNumber] += lambda[iNumberVertice-2 + i*(iNumberVertice-1)+j] * vvCurrentEMatrix[j][k];
                        if (j>0)
                            values[iEleNumber] += lambda[iNumberVertice-2 + iNumberPath*(iNumberVertice-1) + i*(iNumberVertice-1)+j-1] * vvCurrentEMatrix[k][j];
                    }else
                    {
                        values[iEleNumber] = obj_factor * (-2) * vvCurrentSMatrix[j][j];
                        if (j<(iNumberVertice-1))
                            values[iEleNumber] += lambda[iNumberVertice-2 + i*(iNumberVertice-1)+j] * 2 * vvCurrentEMatrix[j][j];
                        if (j>0)
                            values[iEleNumber] += lambda[iNumberVertice-2 + iNumberPath*(iNumberVertice-1) + i*(iNumberVertice-1)+j-1] * 2 * vvCurrentEMatrix[j][j];
                    }
                    iEleNumber++;
                }
        cout<<"eval_h_2: "<<iEleNumber<<endl;
    }*/
    return true;
/*    
  assert (n==4);
  assert (m==3);
  assert(nele_hess==2);
  if(values==NULL)
  {
    iRow[0] = 2;
    jCol[0] = 2;
    
    iRow[1] = 3;
    jCol[1] = 3;
  }
  else {
    values[0] = 2*lambda[0];
    values[1] = 2*lambda[0];
  }
  return true;
  */
}

void
MyTMINLP::finalize_solution(TMINLP::SolverReturn status,
                            Index n, const Number* x, Number obj_value)
{
  std::cout<<"Problem status: "<<status<<std::endl;
  std::cout<<"Objective value: "<<obj_value<<std::endl;
//  if(printSol_ && x != NULL){
    std::cout<<"Solution:"<<std::endl;
    for(int i = 0 ; i < n ; i++){
      std::cout<<"x["<<i<<"] = "<<x[i];
      if(i < n-1) std::cout<<", ";}
    std::cout<<std::endl;
//  }
}
